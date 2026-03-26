/*
 * Copyright Xiangyi Zhang 2021
 * The code may be used for academic, non-commercial purposes only.
 * Please contact me at xiangyi.zhang@polymtl.ca for questions
 * If you have improvements, please contact me!
 */
#include "spp.h"

#include "gurobi_c++.h"

#include <cassert>

#include <algorithm>
#include <list>
#include <set>

/*
The paper section 2.1 (1) and (2)
t_items: all the items
t_idx: the item excluded (idx in the t_items)
flag = false, --> height
flag = true, -->  width
*/

std::ostringstream StripPacking::item::ss;
StripPacking::item::item(const int t_idx, const int t_width, const int t_height)
    : idx(t_idx), width(t_width), height(t_height) {
  subItems.clear();
}

StripPacking::item::item(const int t_idx, const int t_width, const int t_height,
                         const int t_idxHelper)
    : idx(t_idx), width(t_width), height(t_height), idxHelper(t_idxHelper) {
  subItems.clear();
}

std::set<int> StripPacking::computeFX(
    const int t_x, const int t_idx,
    const std::vector<const StripPacking::item*>& t_items, bool flag) {
  std::vector<std::vector<int>> res;
  for (size_t i = 0; i < t_items.size(); ++i) {
    std::vector<int> tmp(t_x + 1, 0);
    res.push_back(tmp);
  }
  for (size_t j = 1; j <= t_x; ++j) res[0][j] = 0;
  for (size_t i = 0; i < t_items.size(); ++i) res[i][0] = 1;
  int itemIdx;
  for (size_t j = 1; j <= t_x; ++j)  // W
  {
    for (size_t i = 1; i < t_items.size(); ++i) {
      if (i - 1 < t_idx)
        itemIdx = i - 1;
      else
        itemIdx = i;
      if (flag) {
        int diff = j - t_items[itemIdx]->width;
        if ((diff) < 0)
          res[i][j] = res[i - 1][j];
        else
          res[i][j] =
              std::max(res[i - 1][j], res[i - 1][j - t_items[itemIdx]->width]);
      } else {
        int diff = j - t_items[itemIdx]->height;
        if ((diff) < 0)
          res[i][j] = res[i - 1][j];
        else
          res[i][j] =
              std::max(res[i - 1][j], res[i - 1][j - t_items[itemIdx]->height]);
      }
    }
  }
  std::set<int> possiblePositions;
  for (size_t j = 0; j <= t_x; ++j) {
    if (res[res.size() - 1][j] == 1) possiblePositions.insert(j);
  }
  return possiblePositions;
}

std::set<int> StripPacking::computeFX(
    const int t_x, const int t_idx,
    const std::vector<StripPacking::item*>& t_items, bool flag) {
  std::vector<std::vector<int>> res;
  for (size_t i = 0; i < t_items.size(); ++i) {
    std::vector<int> tmp(t_x + 1, 0);
    res.push_back(tmp);
  }
  for (size_t j = 1; j <= t_x; ++j) res[0][j] = 0;
  for (size_t i = 0; i < t_items.size(); ++i) res[i][0] = 1;
  int itemIdx;
  for (size_t j = 1; j <= t_x; ++j)  // W
  {
    for (size_t i = 1; i < t_items.size(); ++i) {
      if (i - 1 < t_idx)
        itemIdx = i - 1;
      else
        itemIdx = i;
      if (flag) {
        int diff = j - t_items[itemIdx]->width;
        if ((diff) < 0)
          res[i][j] = res[i - 1][j];
        else
          res[i][j] =
              std::max(res[i - 1][j], res[i - 1][j - t_items[itemIdx]->width]);
      } else {
        int diff = j - t_items[itemIdx]->height;
        if ((diff) < 0)
          res[i][j] = res[i - 1][j];
        else
          res[i][j] =
              std::max(res[i - 1][j], res[i - 1][j - t_items[itemIdx]->height]);
      }
    }
  }
  std::set<int> possiblePositions;
  for (size_t j = 0; j <= t_x; ++j) {
    if (res[res.size() - 1][j] == 1) possiblePositions.insert(j);
  }
  return possiblePositions;
}

int StripPacking::getMaximalHeight(
    const std::vector<const StripPacking::item*>& t_items) {
  int res = -1;
  for (size_t i = 0; i < t_items.size(); ++i) {
    res = std::max(t_items[i]->height, res);
  }
  return res;
}

/*
Build the contiguity parallel machine scheduling problem as a lower bound for
the spp
*/
double StripPacking::solve(
    const std::vector<const StripPacking::item*>& t_allItems,
    const std::map<int, std::set<int>>& t_mapPosWidth,
    const std::map<int, std::set<int>>& t_mapPosHeight, const bool t_Integer) {
  // data preparation
  std::set<int> allPositions;
  for (const auto& it : t_mapPosWidth)
    for (const auto& it2 : it.second) allPositions.insert(it2);
  GRBEnv env(true);
  env.set(GRB_IntParam_OutputFlag, 0);
  env.start();
  GRBModel model(env);
  std::map<std::string, GRBVar> allVars;
  // first constraints set
  for (const auto& it : t_allItems) {
    GRBLinExpr expr;
    for (const auto& it2 : t_mapPosWidth.find(it->idx)->second) {
      auto varName = StripPacking::getVarName(it->idx, it2);
      char varType = t_Integer ? GRB_INTEGER : GRB_CONTINUOUS;
      GRBVar var = model.addVar(0.0, 1.0, 0.0, varType, varName);
      allVars.insert(std::pair<std::string, GRBVar>(varName, var));
      expr += var;
    }
    model.addConstr(expr == 1);
  }
  // second constraints set
  GRBVar z = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "ObjZ");
  for (const auto q : allPositions) {
    GRBLinExpr expr;
    for (const auto it : t_allItems) {
      // calculate W(j, q)
      for (const auto& it2 : t_mapPosWidth.find(it->idx)->second) {
        if (it2 <= q && it2 >= q - it->width + 1) {
          auto iter = allVars.find(StripPacking::getVarName(it->idx, it2));
          assert(iter != allVars.end());
          expr += iter->second * it->height;
        }
      }
    }
    model.addConstr(expr <= z);
  }
  model.setObjective(GRBLinExpr(z), GRB_MINIMIZE);
  // model.write("lowerBound5.lp");

  model.set(GRB_IntParam_Presolve, 2);   // 2=aggressive (closest to CPLEX RepeatPresolve=3 + Reduce=3)
  model.set(GRB_IntParam_Probe, 3);      // 3=aggressive probing, MIP only (mirrors CPLEX MIP::Strategy::Probe=3)
  model.set(GRB_IntParam_Symmetry, 2);   // 2=aggressive, Gurobi max (closest to CPLEX Symmetry=5)
  model.optimize();
  double result = model.get(GRB_DoubleAttr_ObjVal);
  return result;
}

int StripPacking::subSetSum(const std::vector<int>& t_v, const int t_limit)
/*
Args:
        Given a vector of integers and a limit. Find the subset of the integers
of which the sum is the largest but not exceeding the limit This is achieved by
a simple dynamic programming algorithm Returns: The best sum.
*/
{
  int** values = new int*[t_v.size() + 1];
  for (int i = 0; i < t_v.size() + 1; ++i) {
    values[i] = new int[t_limit + 1];
    if (i == 0) {
      for (int j = 0; j < t_limit + 1; ++j) values[i][j] = 0;
    }
    values[i][0] = 0;
  }

  for (int i = 1; i < t_v.size() + 1; ++i) {
    for (int j = 1; j < t_limit + 1; ++j) {
      if (t_v[i - 1] > j)
        values[i][j] = values[i - 1][j];
      else
        values[i][j] = std::max(values[i - 1][j - t_v[i - 1]] + t_v[i - 1],
                                values[i - 1][j]);
    }
  }
  int result = values[t_v.size()][t_limit];
  for (int i = 0; i < t_v.size(); ++i) delete values[i];
  delete values;
  return result;
}

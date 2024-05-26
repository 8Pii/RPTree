#ifndef RPTREE_HPP
#define RPTREE_HPP

#include <Eigen/Dense>
#include <functional>
#include <vector>
#include <memory>
#include <string>
#include <map>

#include "utils.hpp"  


struct TreeNode {
    std::string label;
    std::vector<Eigen::VectorXd> points;
    std::shared_ptr<TreeNode> left;
    std::shared_ptr<TreeNode> right;

    TreeNode(const std::string& lbl, const std::vector<Eigen::VectorXd>& pts, 
             std::shared_ptr<TreeNode> lft, std::shared_ptr<TreeNode> rgt)
        : label(lbl), points(pts), left(lft), right(rgt) {}
};



class Algo_RPTree {
private:
    std::vector<Eigen::VectorXd> data;
    int minSize;
    double c;
    std::string angular;
    double err;
    std::map<std::string, std::pair<Eigen::VectorXd, std::vector<Eigen::VectorXd>>> dic;

    std::shared_ptr<TreeNode> MakeTree(std::vector<Eigen::VectorXd>& S, int depth = 0, const std::string& label = "root");
    std::function<bool(const Eigen::VectorXd&)> ChooseRule(const std::vector<Eigen::VectorXd>& S);
    std::pair<Eigen::VectorXd, double> compute_centroid(const std::vector<Eigen::VectorXd>& pts);
    void process_tree(const std::shared_ptr<TreeNode>& node);
    int countLeafNodes(const std::shared_ptr<TreeNode>& node);

public:
    Algo_RPTree(const std::vector<Eigen::VectorXd>& data, int minSize, double c, const std::string& angular);
    void compute_error();  // Public interface to start error computation
    double get_error() const;
};

#endif // RPTREE_HPP
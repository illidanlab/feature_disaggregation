# KDD2018_FeatureDisaggregation
## Summary
This is the repo for the source code of KDD paper: [Enhancing Predictive Modeling of Nested Spatial Data through
Group-Level Feature Disaggregation](http://delivery.acm.org/10.1145/3230000/3220091/p1784-liu.pdf?ip=69.63.237.24&id=3220091&acc=OPENTOC&key=4D4702B0C3E38B35%2E4D4702B0C3E38B35%2E4D4702B0C3E38B35%2E054E54E275136550&__acm__=1538011558_ecbac57817a93ca74977e4a33f16231c).

Feature disaggregation is designed for dealing with nested spatial data. When the data can be divided to different task and there are nested structure among the features, feature disaggregation can give high-resolution group-level feature by making two assumptions: 1). the disaggregated values should keep spatial contiguity, 2). the disaggregated values should not be too far away from original group-level feature.

## Compatibility
This code is based on MATLAB_R2016a. Also, it is necessary to download [cvx](http://cvxr.com/cvx/) and [MALSAR](https://github.com/jiayuzhou/MALSAR) to run the baselines.

## Reference
Liu, B., Tan, P. N., & Zhou, J. (2018, July). Enhancing Predictive Modeling of Nested Spatial Data through Group-Level Feature Disaggregation. In Proceedings of the 24th ACM SIGKDD International Conference on Knowledge Discovery & Data Mining (pp. 1784-1793). ACM.

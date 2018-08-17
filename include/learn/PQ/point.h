// Basic storage class for an n-dimensional point with a cluster assignment.
//
// Author: Felix Duvallet

#ifndef __KMEANS_POINT_H__
#define __KMEANS_POINT_H__

#include <vector>
#include <iostream>

class Point {
 public:
  Point() { }

  // Initialize the number of dimensions, optionally set all values to zero.
  Point(int num_dimensions, bool init_zeros = true);

  Point(float x, float y, float z);

  // Initialize from a vector.
  Point(const std::vector<float> &vector);

  ~Point() { }

  // Compute distance between two points.
  static float distance(const Point &p1, const Point &p2);

  // Adds a point to the current point.
  void add(const Point &point);

  // Update the cluster assignment. Returns true if cluster assignment
  // changed, false if it stayed the same.
  bool update(int k);

  // Members: the data, the number of dimensions, and the cluster ID.
  std::vector<float> data_;
  int dimensions_;
  int cluster_;

  friend std::ostream &operator<<(std::ostream &target, const Point &point);
};
Point::Point(int num_dimensions, bool init_zeros)
  : cluster_(-1),
    dimensions_(num_dimensions) {
  if (init_zeros) {  // default is true.
    for (int idx = 0; idx < dimensions_; ++idx) {
      data_.push_back(0.0);
    }
  }
}

/*Point::Point(float x, float y, float z)
  : Point(3, false) {
  data_.clear();
  data_.push_back(x);
  data_.push_back(y);
  data_.push_back(z);
}*/

Point::Point(const std::vector<float> &vector)
  : cluster_(-1) {
  dimensions_ = (int) vector.size();
  data_.clear();
  data_.insert(data_.begin(), vector.begin(), vector.end());
}

bool Point::update(int k) {
  const bool ret = cluster_ != k;
  cluster_ = k;
  return ret;
}

// static
float Point::distance(const Point &p1, const Point &p2) {
  assert(p1.dimensions_ == p2.dimensions_);
  float dist = 0.0;

  for (int idx = 0; idx < p1.dimensions_; ++idx) {
    const float tmp = p1.data_[idx] - p2.data_[idx];
    dist += tmp * tmp;
  }
  return sqrt(dist);
}

void Point::add(const Point &point) {
  assert(dimensions_ == point.dimensions_);
  for (int idx = 0; idx < dimensions_; ++idx) {
    data_[idx] += point.data_[idx];
  }
}

std::ostream &operator<<(std::ostream &target, const Point &point) {
  target << "[";
  for (const float &d : point.data_) {
    target << d << ", ";
  }
  target << "]";
  return target;
}


#endif  // __KMEANS_POINT_H_

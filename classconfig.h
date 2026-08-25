#pragma once

namespace cc {

struct cell_class;
struct face_class;
struct node_class;

struct node_class{
    int number = 0;
    double x = 0.0;
    double y = 0.0;

    node_class() = default;
    
    // 记录点的坐标
    node_class(int number_,double x_,double y_);
};



}
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "classconfig.h"

static bool ifin_node(int num[], int number) {
  for (int i = 0; i < 4; i++) {if (number == num[i]) {return true;}}return false;
}

void findnode(cc::cell_class &cell) {
  short place = 0;
  for (int i = 0; i < cell.ecnt; i++) {
    if (!ifin_node(cell.node, cell.nei[i]->node.x - 1)) {
      cell.node[place] = cell.nei[i]->node.x - 1;
      place++;
    }
    if (!ifin_node(cell.node, cell.nei[i]->node.y - 1)) {
      cell.node[place] = cell.nei[i]->node.y - 1;
      place++;
    }
    if(place == cell.ecnt){break;}
  }
}

static double triangle_area(double x1, double y1, double x2, double y2, double x3, double y3) {
    return std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) * 0.5;
}

static cc::vec2 triangle_center(double x1, double y1, double x2, double y2, double x3, double y3){
    return {(x1+x2+x3)/3,(y1+y2+y3)/3};
}

double volume(cc::cell_class &cell){
    const auto& p0 = cc::NodeList[cell.node[0]];
    const auto& p1 = cc::NodeList[cell.node[1]];
    const auto& p2 = cc::NodeList[cell.node[2]];
    if(cell.ecnt == 3){
        cell.vol = triangle_area(p0.x,p0.y,p1.x,p1.y,p2.x,p2.y);
        return cell.vol;
    }
    const auto& p3 = cc::NodeList[cell.node[3]];
    double S[4];
    S[0] = triangle_area(p0.x,p0.y,p1.x,p1.y,p2.x,p2.y);
    S[1] = triangle_area(p0.x,p0.y,p1.x,p1.y,p3.x,p3.y);
    S[2] = triangle_area(p0.x,p0.y,p2.x,p2.y,p3.x,p3.y);
    S[3] = triangle_area(p1.x,p1.y,p2.x,p2.y,p3.x,p3.y);
    cell.vol = (S[0]+S[1]+S[2]+S[3])*0.5;
    return cell.vol;
}

void center(cc::cell_class &cell){
    const auto& p0 = cc::NodeList[cell.node[0]];
    const auto& p1 = cc::NodeList[cell.node[1]];
    const auto& p2 = cc::NodeList[cell.node[2]];
    if(cell.ecnt == 3){
        cell.center = triangle_center(p0.x,p0.y,p1.x,p1.y,p2.x,p2.y);
        return;
    }
    const auto& p3 = cc::NodeList[cell.node[3]];
    double S[4];cc::vec2 G[4];

    S[0] = triangle_area(p0.x,p0.y,p1.x,p1.y,p2.x,p2.y);
    G[0] = triangle_center(p0.x,p0.y,p1.x,p1.y,p2.x,p2.y);
    S[1] = triangle_area(p0.x,p0.y,p1.x,p1.y,p3.x,p3.y);
    G[1] = triangle_center(p0.x,p0.y,p1.x,p1.y,p3.x,p3.y);
    S[2] = triangle_area(p0.x,p0.y,p2.x,p2.y,p3.x,p3.y);
    G[2] = triangle_center(p0.x,p0.y,p2.x,p2.y,p3.x,p3.y);
    S[3] = triangle_area(p1.x,p1.y,p2.x,p2.y,p3.x,p3.y);
    G[3] = triangle_center(p1.x,p1.y,p2.x,p2.y,p3.x,p3.y);

    cell.center = {(S[0]*G[0].x  + S[1]*G[1].x  + S[2]*G[2].x  + S[3]*G[3].x)/(S[0]+S[1]+S[2]+S[3]),
                    (S[0]*G[0].y + S[1]*G[1].y + S[2]*G[2].y + S[3]*G[3].y)/(S[0]+S[1]+S[2]+S[3])};
    
}

void walldistance_alternative(cc::cell_class &cell,double H){
  cell.sad = cell.center.y < (H - cell.center.y) ? cell.center.y : H - cell.center.y;
}
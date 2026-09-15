#include "readmesh.h"
#include "classconfig.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>

namespace {

std::ifstream input;
int line_no = 0;

bool mesh_fail(const std::string& msg){
    fprintf(stderr,"Error: mesh line %d: %s\n",line_no,msg.c_str());
    return false;
}

bool next_line(std::string& text){
    if(!std::getline(input,text)){
        return mesh_fail("unexpected end of mesh");
    }
    line_no++;
    if(!text.empty() && text.back() == '\r'){
        text.pop_back();
    }
    return true;
}

bool expect(const char* wanted){
    std::string text;
    if(!next_line(text)){
        return false;
    }
    if(text != wanted){
        return mesh_fail(std::string("expected ") + wanted);
    }
    return true;
}

template <typename... value_type>
bool parse_row(const std::string& text,value_type&... values){
    std::istringstream row(text);
    if(!(row >> ... >> values)){
        return mesh_fail("invalid mesh row");
    }
    row >> std::ws;
    if(!row.eof()){
        return mesh_fail("unexpected data after mesh row");
    }
    return true;
}

}

bool readmesh(const char* path){
    input.open(path);
    if(!input){
        fprintf(stderr,"Error: Cannot open mesh: %s\n",path);
        return false;
    }
    std::string line;
    int node_count = 0,type_count = 0;
    if(!next_line(line)) return false;
    if(!parse_row(line,node_count,cc::face_num,cc::cell_num,type_count)) return false;
    if(node_count < 4 || cc::face_num < 4 || cc::cell_num < 1 || type_count != 3){
        return mesh_fail("expected positive counts and exactly INTER/WALL/FAR groups");
    }
    const std::map<std::string,short> supported = {
        {"INTER",cc::INTER},{"WALL",cc::WALL},{"FAR",cc::FAR}};
    std::vector<std::pair<std::string,short>> groups;
    std::set<std::string> names,kinds;
    for(int i=0;i<type_count;i++){
        if(!next_line(line)) return false;
        std::string name,equal,kind;
        if(!parse_row(line,name,equal,kind)) return false;
        if(equal != "=" || !supported.count(kind) || !names.insert(name).second ||
           !kinds.insert(kind).second){
            return mesh_fail("invalid or duplicate boundary group");
        }
        groups.emplace_back(name,supported.at(kind));
    }
    if(!expect("(node)")) return false;
    cc::NodeList.reserve(node_count);
    for(int i=1;i<=node_count;i++){
        if(!next_line(line)) return false;
        int index = 0;
        double x = 0.0,y = 0.0;
        if(!parse_row(line,index,x,y)) return false;
        if(index != i || !std::isfinite(x) || !std::isfinite(y)){
            return mesh_fail("invalid node index or coordinates");
        }
        cc::NodeList.emplace_back(index,x,y);
    }
    if(!expect("(end)")) return false;
    if(!expect("(edge)")) return false;
    cc::FaceList.reserve(cc::face_num);
    std::set<int> face_ids;
    for(const auto& group : groups){
        if(!expect(group.first.c_str())) return false;
        for(;;){
            if(!next_line(line)) return false;
            if(line == "(end)"){
                break;
            }
            int index = 0,a = 0,b = 0,left = 0,right = 0;
            if(!parse_row(line,index,a,b,left,right)) return false;
            if(index < 1 || index > cc::face_num || !face_ids.insert(index).second){
                return mesh_fail("invalid or duplicate face index");
            }
            if(a < 1 || a > node_count || b < 1 || b > node_count || a == b){
                return mesh_fail("invalid face nodes");
            }
            if(left < 0 || left > cc::cell_num || right < 0 || right > cc::cell_num ||
               left == right){
                return mesh_fail("invalid adjacent cell indices");
            }
            const bool interior = left != 0 && right != 0;
            if(interior != (group.second == cc::INTER)){
                return mesh_fail("boundary type disagrees with adjacency");
            }
            const cc::node_class& pa = cc::NodeList[a-1];
            const cc::node_class& pb = cc::NodeList[b-1];
            if(!(std::hypot(pa.x - pb.x,pa.y - pb.y) > 0.0)){
                return mesh_fail("zero-length face");
            }
            cc::FaceList.emplace_back(index,a,b,left,right,group.second);
        }
    }
    if(cc::FaceList.size() != static_cast<size_t>(cc::face_num)){
        return mesh_fail("face count mismatch");
    }
    std::sort(cc::FaceList.begin(),cc::FaceList.end(),
              [](const cc::face_class& a,const cc::face_class& b){ return a.index < b.index; });
    if(!expect("(end)")) return false;
    if(!expect("(cell)")) return false;
    cc::CellList.reserve(cc::cell_num);
    std::vector<int> uses(cc::face_num,0);
    for(int i=1;i<=cc::cell_num;i++){
        if(!next_line(line)) return false;
        int index = 0;
        std::array<int,4> faces{};
        if(!parse_row(line,index,faces[0],faces[1],faces[2],faces[3])) return false;
        if(index != i || std::set<int>(faces.begin(),faces.end()).size() != 4){
            return mesh_fail("expected sequential quadrilateral cells with four distinct faces");
        }
        std::map<int,int> degree;
        for(int id : faces){
            if(id < 1 || id > cc::face_num){
                return mesh_fail("cell face index out of range");
            }
            const cc::face_class& face = cc::FaceList[id-1];
            if(face.cell_1 != i && face.cell_2 != i){
                return mesh_fail("cell/face adjacency mismatch");
            }
            ++degree[face.node[0]->number];
            ++degree[face.node[1]->number];
            ++uses[id-1];
        }
        if(degree.size() != 4 ||
           std::any_of(degree.begin(),degree.end(),
                       [](const std::pair<const int,int>& entry){ return entry.second != 2; })){
            return mesh_fail("cell edges do not form a quadrilateral");
        }
        for(int id : faces){
            const cc::face_class& edge = cc::FaceList[id-1];
            const cc::node_class& a = *edge.node[0];
            const cc::node_class& b = *edge.node[1];
            double side = 0.0;
            for(const auto& entry : degree){
                if(entry.first == a.number || entry.first == b.number){
                    continue;
                }
                const cc::node_class& c = cc::NodeList[entry.first-1];
                const double cross = (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x);
                if(!std::isfinite(cross) || cross == 0.0 ||
                   (side != 0.0 && std::signbit(side) != std::signbit(cross))){
                    return mesh_fail("cell is degenerate, concave or self-intersecting");
                }
                side = cross;
            }
        }
        cc::CellList.emplace_back(index,faces[0],faces[1],faces[2],faces[3]);
    }
    for(const cc::face_class& face : cc::FaceList){
        const int wanted = face.type == cc::INTER ? 2 : 1;
        if(uses[face.index-1] != wanted){
            return mesh_fail("face reference count mismatch");
        }
    }
    input >> std::ws;
    if(!input.eof()){
        return mesh_fail("unexpected data after cell section");
    }
    printf("Mesh: nodes=%d faces=%d cells=%d\n",node_count,cc::face_num,cc::cell_num);
    return true;
}

bool linkmesh(){
    for(cc::cell_class& cell : cc::CellList){
        for(int i=0;i<cell.ecnt;i++){
            cell.faces[i] = cc::link_face(cell.face[i]);
        }
    }
    for(cc::face_class& face : cc::FaceList){
        face.nei[0] = cc::link_cell(face.cell_1);
        face.nei[1] = cc::link_cell(face.cell_2);
        if(face.type == cc::WALL){
            cc::WallFaces.push_back(&face);
        }
        if(face.type == cc::FAR){
            cc::FarFaces.push_back(&face);
        }
    }
    for(cc::cell_class& cell : cc::CellList){
        for(int i=0;i<cell.ecnt;i++){
            const cc::face_class& face = *cell.faces[i];
            cell.nei[i] = face.nei[0] == &cell ? face.nei[1] : face.nei[0];
        }
    }
    if(cc::WallFaces.empty() || cc::FarFaces.empty()){
        fprintf(stderr,"Error: mesh needs nonempty wall and farfield boundaries\n");
        return false;
    }
    printf("Boundaries: wall=%zu farfield=%zu\n",cc::WallFaces.size(),cc::FarFaces.size());
    return true;
}

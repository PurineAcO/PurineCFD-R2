#include "config.h"
#include "timarch.h"
#include "physic.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>

namespace {

// 极简JSON解析器, 支持//与/* */注释
struct JVal {
    enum Type{NUL,BOL,NUM,STR,OBJ,ARR} t = NUL;
    bool b = false; double n = 0.0; std::string s;
    std::map<std::string,JVal> o; std::vector<JVal> a;
    const JVal* find(const std::string& k) const {
        auto it = o.find(k);
        return (it == o.end()) ? nullptr : &it->second;
    }
};

class Parser {
    const std::string& s; std::size_t i = 0;
    void ws(){
        while(i < s.size()){
            char c = s[i];
            if(c==' '||c=='\t'||c=='\n'||c=='\r'){ i++; }
            else if(c=='/' && i+1<s.size() && s[i+1]=='/'){ while(i<s.size()&&s[i]!='\n') i++; }
            else if(c=='/' && i+1<s.size() && s[i+1]=='*'){ i+=2; while(i+1<s.size() && !(s[i]=='*'&&s[i+1]=='/')) i++; i+=2; }
            else break;
        }
    }
    std::string str(){
        i++; std::string out;
        while(i < s.size()){
            char c = s[i++];
            if(c=='"') break;
            if(c=='\\' && i<s.size()){
                char e = s[i++];
                switch(e){
                    case 'n': out+='\n'; break; case 't': out+='\t'; break;
                    case 'r': out+='\r'; break; case '"': out+='"'; break;
                    case '\\': out+='\\'; break; case '/': out+='/'; break;
                    case 'b': out+='\b'; break; case 'f': out+='\f'; break;
                    default: out+=e;
                }
            }else out+=c;
        }
        return out;
    }
    JVal num(){
        std::size_t st = i;
        while(i<s.size()){
            char c = s[i];
            if(c=='-'||c=='+'||c=='.'||c=='e'||c=='E'||(c>='0'&&c<='9')) i++; else break;
        }
        JVal v; v.t=JVal::NUM; v.n = std::atof(s.substr(st,i-st).c_str());
        return v;
    }
    JVal obj(){
        i++; JVal v; v.t=JVal::OBJ;
        while(true){
            ws(); if(i>=s.size()) break;
            if(s[i]=='}'){ i++; break; }
            if(s[i] != '"'){ i++; continue; }
            std::string k = str();
            ws(); if(i<s.size() && s[i]==':') i++;
            v.o[k] = value();
            ws(); if(i<s.size() && s[i]==',') i++;
        }
        return v;
    }
    JVal arr(){
        i++; JVal v; v.t=JVal::ARR;
        while(true){
            ws(); if(i>=s.size()) break;
            if(s[i]==']'){ i++; break; }
            v.a.push_back(value());
            ws(); if(i<s.size() && s[i]==',') i++;
        }
        return v;
    }
public:
    explicit Parser(const std::string& src):s(src){}
    JVal value(){
        ws(); if(i>=s.size()) return JVal();
        char c = s[i];
        if(c=='{') return obj();
        if(c=='[') return arr();
        if(c=='"'){ JVal v; v.t=JVal::STR; v.s=str(); return v; }
        if(c=='t'){ i+=4; JVal v; v.t=JVal::BOL; v.b=true; return v; }
        if(c=='f'){ i+=5; JVal v; v.t=JVal::BOL; v.b=false; return v; }
        if(c=='n'){ i+=4; return JVal(); }
        return num();
    }
};

double dbl(const JVal* v, const char* k, double def){
    if(!v) return def;
    const JVal* c = v->find(k);
    return (c && c->t==JVal::NUM) ? c->n : def;
}
bool flg(const JVal* v, const char* k, bool def){
    if(!v) return def;
    const JVal* c = v->find(k);
    return (c && c->t==JVal::BOL) ? c->b : def;
}

} // namespace

void config::load(const char* path){
    std::ifstream f(path);
    if(!f.is_open()){ std::printf("[config] 无法打开 %s\n", path); return; }
    std::stringstream ss; ss << f.rdbuf();
    std::string text = ss.str();
    Parser ps(text);
    JVal root = ps.value();

    const JVal* io = root.find("io");
    static std::string sMesh, sLog, sField;
    if(io){
        const JVal* m = io->find("mesh"); if(m && m->t==JVal::STR){ sMesh=m->s; cc::meshpath=sMesh.c_str(); }
        const JVal* g = io->find("log");  if(g && g->t==JVal::STR){ sLog=g->s;  cc::testpath=sLog.c_str(); }
        const JVal* d = io->find("field");if(d && d->t==JVal::STR){ sField=d->s;cc::fieldpath=sField.c_str(); }
    }

    const JVal* tm = root.find("time");
    if(tm){
        cc::max_step    = (long long)dbl(tm,"max_step",(double)cc::max_step);
        cc::total_time  = dbl(tm,"total_time",cc::total_time);
        fatime::CFL     = dbl(tm,"cfl",fatime::CFL);
        fatime::USE_GLOBAL_DT = flg(tm,"global_dt",fatime::USE_GLOBAL_DT);
        config::dump_step = (int)dbl(tm,"dump_step",(double)config::dump_step);
        config::conv_step = (int)dbl(tm,"conv_step",(double)config::conv_step);
        const JVal* rk = tm->find("rk_coeff");
        if(rk && rk->t==JVal::ARR && !rk->a.empty()){
            RK::RK.clear();
            for(const JVal& c : rk->a) if(c.t==JVal::NUM) RK::RK.push_back(c.n);
        }
    }

    const JVal* bc = root.find("boundary");
    if(bc){
        const JVal* inlet = bc->find("inlet");
        if(inlet){
            cc::VIL_DEFINE.u = dbl(inlet,"u",cc::VIL_DEFINE.u);
            cc::VIL_DEFINE.v = dbl(inlet,"v",cc::VIL_DEFINE.v);
            cc::VIL_DEFINE.T = dbl(inlet,"T",cc::VIL_DEFINE.T);
            cc::VIL_DEFINE.p = dbl(inlet,"p",cc::VIL_DEFINE.p);
        }
        const JVal* far = bc->find("farfield");
        if(far){
            double p = dbl(far,"p",101325.0);
            double T = dbl(far,"T",300.0);
            double Ma= dbl(far,"Ma",0.0);
            double aoa= dbl(far,"AOA_deg",0.0);
            double U  = Ma * get_sonic_velocity(T);
            double u  = U * cos(deg2rad(aoa));
            double v  = U * sin(deg2rad(aoa));
            cc::FAR_DEFINE = {u,v,T,p};
        }
        const JVal* out = bc->find("outlet");
        if(out){
            cc::POL_DEFINE.p = dbl(out,"p",cc::POL_DEFINE.p);
            cc::POL_DEFINE.T = dbl(out,"T",cc::POL_DEFINE.T);
        }
    }

    // 湍流/粘性
    const JVal* vis = root.find("viscous");
    static std::string sTurb;
    if(vis){
        cc::viscous = flg(vis,"ifviscous",false);
        const JVal* mod = vis->find("model");
        if(mod && mod->t==JVal::STR){ sTurb = mod->s; cc::turb_model = sTurb.c_str(); }
    }

    std::printf("[config] mesh=%s CFL=%.2f max_step=%lld viscous=%s model=%s\n",
        cc::meshpath, fatime::CFL, cc::max_step,
        cc::viscous ? "ON" : "OFF", cc::turb_model ? cc::turb_model : "-");
}

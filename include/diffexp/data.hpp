#pragma once
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace diffexp::data {
// Data-only reader for published ancillary expressions. No evaluator, kernel,
// shell expansion, assignments or executable Wolfram constructs are supported.
struct Expr {
  std::string head;
  std::vector<Expr> args;
  bool atom() const { return args.empty(); }
};

class Reader {
 public:
  explicit Reader(std::string input):text_(std::move(input)) {next();}
  Expr read() {
    auto e=rule(); if(!token_.empty()) fail("trailing data"); return e;
  }
 private:
  std::string text_,token_; std::size_t pos_=0;
  [[noreturn]] void fail(const char* why) const {
    throw std::invalid_argument(std::string("ancillary data: ")+why+" at byte "+std::to_string(pos_));
  }
  void next() {
    while(pos_<text_.size()) {
      if(std::isspace(static_cast<unsigned char>(text_[pos_]))) {++pos_;continue;}
      if(text_.compare(pos_,2,"\\\n")==0) {pos_+=2;continue;}
      if(text_.compare(pos_,2,"(*")==0) {
        pos_+=2; unsigned depth=1;
        while(depth && pos_<text_.size()) {
          if(text_.compare(pos_,2,"(*")==0){++depth;pos_+=2;}
          else if(text_.compare(pos_,2,"*)")==0){--depth;pos_+=2;}
          else ++pos_;
        }
        if(depth) fail("unterminated comment");
        continue;
      }
      break;
    }
    if(pos_==text_.size()) {token_.clear();return;}
    const auto start=pos_; const char c=text_[pos_++];
    if(std::isdigit(static_cast<unsigned char>(c)) || c=='.') {
      while(pos_<text_.size()) {
        char t=text_[pos_];
        if(std::isdigit(static_cast<unsigned char>(t)) || t=='.' || t=='`') ++pos_;
        else if(text_.compare(pos_,2,"\\\n")==0) pos_+=2;
        else if(text_.compare(pos_,2,"*^")==0) {
          pos_+=2; if(pos_<text_.size() && (text_[pos_]=='-' || text_[pos_]=='+')) ++pos_;
        } else break;
      }
    } else if(std::isalpha(static_cast<unsigned char>(c)) || c=='$' || c=='_') {
      while(pos_<text_.size() && (std::isalnum(static_cast<unsigned char>(text_[pos_])) || text_[pos_]=='_' || text_[pos_]=='$' || text_[pos_]=='`')) ++pos_;
    } else if(c=='-' && pos_<text_.size() && text_[pos_]=='>') ++pos_;
    token_=text_.substr(start,pos_-start);
    for(std::size_t p;(p=token_.find("\\\n"))!=std::string::npos;) token_.erase(p,2);
  }
  bool take(const std::string& token) {if(token_!=token)return false;next();return true;}
  void need(const std::string& token) {if(!take(token))fail("unexpected token");}
  Expr rule() {
    auto e=add(); if(take("->")) return {"Rule",{std::move(e),rule()}}; return e;
  }
  Expr add() {
    auto e=product();
    while(token_=="+" || token_=="-") {auto op=token_;next();e={op,{std::move(e),product()}};}
    return e;
  }
  bool primary_start() const {
    if(token_.empty())return false;
    return std::isalnum(static_cast<unsigned char>(token_[0])) || token_[0]=='.' || token_=="(";
  }
  Expr product() {
    auto e=unary();
    while(token_=="*" || token_=="/" || primary_start()) {
      auto op=token_=="/"?"/":"*";
      if(token_=="/" || token_=="*") next();
      e={op,{std::move(e),unary()}};
    }
    return e;
  }
  Expr unary() {
    if(take("+"))return unary();
    if(take("-"))return {"neg",{unary()}};
    auto e=primary();
    if(take("^"))return {"^",{std::move(e),unary()}};
    return e;
  }
  Expr primary() {
    if(take("(")) {auto e=rule();need(")");return e;}
    if(take("{")) return list("List","}");
    if(token_.empty() || (!std::isalnum(static_cast<unsigned char>(token_[0])) && token_[0]!='.' && token_[0]!='$' && token_[0]!='_'))fail("expected expression");
    auto head=token_;next();
    if(take("["))return list(head,"]");
    return {head,{}};
  }
  Expr list(std::string head,const std::string& close) {
    Expr out{std::move(head),{}};
    if(take(close))return out;
    do {out.args.push_back(rule());}while(take(","));
    need(close);return out;
  }
};

inline Expr read_file(const std::string& path) {
  std::ifstream in(path); if(!in)throw std::runtime_error("cannot read ancillary file: "+path);
  std::ostringstream text;text<<in.rdbuf();return Reader(text.str()).read();
}
inline bool number(const Expr& e) {
  return e.atom() && !e.head.empty() && (std::isdigit(static_cast<unsigned char>(e.head[0])) || e.head[0]=='.');
}
inline long integer(const Expr& e) {
  if(e.head=="neg" && e.args.size()==1)return -integer(e.args[0]);
  if(!number(e))throw std::invalid_argument("expected literal integer");
  std::size_t used=0;auto out=std::stol(e.head,&used);
  if(used!=e.head.size())throw std::invalid_argument("expected literal integer");
  return out;
}
}  // namespace diffexp::data

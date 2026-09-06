#include "diffexp/ft_spectral_checkpoint.hpp"
#include <iostream>
using namespace diffexp;using B=Jet::Ball;namespace cp=ft_spectral_checkpoint;
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
template<class F>void rejects(F f,const char* message){try{f();}catch(const std::exception&){return;}throw std::runtime_error(message);}
void equal(const LaurentRows& a,const LaurentRows& b){require(a.low==b.low&&a.high==b.high&&a.columns()==b.columns()&&a.coefficients.size()==b.coefficients.size(),"checkpoint shape round trip");for(unsigned r=0;r<a.coefficients.size();++r)for(unsigned i=0;i<a.columns();++i)for(unsigned e=0;e<a.coefficients[r][i].size();++e)require(acb_equal(a.coefficients[r][i][e].raw(),b.coefficients[r][i][e].raw()),"checkpoint changed exact ball bits");}
int main(){try {
 B::set_precision(256);ExactField field({"x","eps","I"});Exact x(field,"x"),eps(field,"eps"),zero(field,0),one(field,1);
 ExactEpsilonMatrix matrix{{zero}},forcing{{one/eps}};LaurentRows initial{0,1,{{{B(1),B(0)}}}};arb_add_error_2exp_si(acb_realref(initial.coefficients[0][0][0].raw()),-200);
 std::vector<Exact> path{zero,one/one.constant(2),one};ft_spectral::Options spectral;spectral.accuracy_goal=20;AdjointOptions native;ft_spectral::Diagnostics diagnostics;
 const auto directory=std::filesystem::temp_directory_path()/("diffexp-ft-spectral-checkpoint-"+std::to_string(::getpid()));
 struct Cleanup{std::filesystem::path p;~Cleanup(){std::error_code error;std::filesystem::remove_all(p,error);}} cleanup{directory};
 std::size_t reused=0;cp::Storage storage{directory,64*1024*1024,&reused};
 auto key=cp::identity(matrix,initial,forcing,path,spectral,native);auto file=directory/(key+".json");
 auto first=cp::try_transport(matrix,initial,forcing,path,spectral,diagnostics,native,storage);require(bool(first),diagnostics.reason.c_str());require(std::filesystem::exists(file)&&reused==0,"successful arm was not published");
 auto second=cp::try_transport(matrix,initial,forcing,path,spectral,diagnostics,native,storage);require(bool(second)&&reused==1,"completed arm not reused");equal(*first,*second);require(diagnostics.legs==2&&diagnostics.factorizations==0&&diagnostics.numerical_seconds==0,"cache reuse repeated or reported numerical work");
 auto envelope=adjoint_checkpoint::detail::read(file,storage.max_bytes);
 B::set_precision(64);ft_spectral::Diagnostics loaded;auto exact=cp::detail::decode(envelope,key,initial,-1,2,spectral,loaded);B::set_precision(256);equal(*first,exact);
 auto changed=[&](auto change){auto options=spectral;change(options);require(cp::identity(matrix,initial,forcing,path,options,native)!=key,"spectral option omitted from identity");};
 changed([](auto& o){o.diagonal_gauge=!o.diagonal_gauge;});changed([](auto& o){o.endpoint_clustering=!o.endpoint_clustering;});changed([](auto& o){o.conservative=!o.conservative;});changed([](auto& o){++o.accuracy_goal;});changed([](auto& o){o.max_nodes=96;});changed([](auto& o){++o.max_block_size;});changed([](auto& o){++o.max_block_nodes;});changed([](auto& o){++o.max_cells;});changed([](auto& o){o.seconds_budget+=1;});
 auto alternate=native;++alternate.taylor_order;require(cp::identity(matrix,initial,forcing,path,spectral,alternate)!=key,"native settings omitted from identity");
 auto boundary=initial;boundary.coefficients[0][0][0]+=B(1);require(cp::identity(matrix,boundary,forcing,path,spectral,native)!=key,"boundary omitted from identity");
 require(cp::identity({{one}},initial,forcing,path,spectral,native)!=key,"matrix omitted from identity");require(cp::identity(matrix,initial,{{one}},path,spectral,native)!=key,"forcing omitted from identity");
 require(cp::identity(matrix,initial,forcing,{zero,one},spectral,native)!=key,"path omitted from identity");B::set_precision(320);require(cp::identity(matrix,initial,forcing,path,spectral,native)!=key,"working bits omitted from identity");B::set_precision(256);
 auto tiny=storage;tiny.max_bytes=1;rejects([&]{cp::try_transport(matrix,initial,forcing,path,spectral,diagnostics,native,tiny);},"file read budget ignored");
 auto publish=[&](boost::json::value payload){auto bytes=artifacts::detail::canonical(boost::json::object{{"payload",payload},{"sha256",artifacts::detail::sha256(artifacts::detail::canonical(payload))}});adjoint_checkpoint::detail::publish(file,bytes);};
 auto payload=envelope.as_object().at("payload");auto bad=payload;bad.as_object()["identity"]="wrong";publish(bad);rejects([&]{cp::try_transport(matrix,initial,forcing,path,spectral,diagnostics,native,storage);},"matching file identity corruption silently restarted");
 bad=payload;bad.as_object().at("rows").as_object()["low"]=0;publish(bad);rejects([&]{cp::try_transport(matrix,initial,forcing,path,spectral,diagnostics,native,storage);},"corrupt matching row window accepted");
 bad=payload;bad.as_object().at("diagnostics").as_object()["legs"]=1;publish(bad);rejects([&]{cp::try_transport(matrix,initial,forcing,path,spectral,diagnostics,native,storage);},"partial arm accepted");
 auto nonfinite=*first;acb_indeterminate(nonfinite.coefficients[0][0][0].raw());bad=payload;bad.as_object()["rows"]=numerical_rows_io::exact_rows(nonfinite);publish(bad);rejects([&]{cp::try_transport(matrix,initial,forcing,path,spectral,diagnostics,native,storage);},"nonfinite checkpoint accepted");
 auto corrupt=envelope;corrupt.as_object()["sha256"]="invalid";adjoint_checkpoint::detail::publish(file,artifacts::detail::canonical(corrupt));rejects([&]{cp::try_transport(matrix,initial,forcing,path,spectral,diagnostics,native,storage);},"checksum corruption silently restarted");
 publish(payload);auto coarse=spectral;coarse.max_nodes=16;coarse.accuracy_goal=40;auto rejection=cp::try_transport({{-one}},initial,{{zero}},path,coarse,diagnostics,native,storage);require(!rejection,"underresolved spectral arm accepted");require(!std::filesystem::exists(directory/(cp::identity({{-one}},initial,{{zero}},path,coarse,native)+".json")),"spectral rejection published a checkpoint");
 // Publication budget failure must not leave a completed matching file.
 auto new_path=std::vector<Exact>{zero,one};rejects([&]{cp::try_transport(matrix,initial,forcing,new_path,spectral,diagnostics,native,tiny);},"publication budget ignored");require(!std::filesystem::exists(directory/(cp::identity(matrix,initial,forcing,new_path,spectral,native)+".json")),"oversized checkpoint published");
 std::cout<<"FT spectral checkpoint identity, exact balls, corruption and completed-arm reuse passed\n";
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}

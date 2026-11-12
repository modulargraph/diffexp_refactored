#include "diffexp/finite_masters.hpp"
#include <iostream>
using namespace diffexp;
namespace fm=finite_masters;
using Matrix=ExactEpsilonMatrix;
void check(bool b,const char* message){if(!b)throw std::runtime_error(message);}
int main(){try {
  ExactField field({"x","eps","I"});
  auto e=[&](const char* s){return Exact(field,s);};
  const Exact x(field,"x");
  // The convergence test on one-loop families with two lines.
  {
    const auto family=feynman::example_family("bubble");
    const ibp::PropagatorBasis massive(ibp::quadratic_family(family.momenta,x,family.physical_count));
    const auto slots=massive.space.size();
    const auto integral=[&](int a,int b){ibp::Integral out(slots,0);out[0]=a;out[1]=b;return out;};
    const fm::Convergence two(massive,1,Rational(2)),four(massive,1,Rational(4));
    check(two.minimal_total()==2 && four.minimal_total()==3,"smallest total with a finite Gamma function");
    check(two(integral(1,1)) && !two(integral(1,0)) && two(integral(2,0)),"massive lines in two dimensions");
    check(!four(integral(1,1)) && four(integral(2,1)) && four(integral(1,2)) && !four(integral(2,0)) && four(integral(3,0)),"massive lines in four dimensions");
    auto numerator=integral(2,1);numerator.back()=slots>2?-1:numerator.back();
    check(slots==2 || !four(numerator),"numerators are not classified");
    const auto found=fm::convergent_integrals(four,integral(1,1),4);
    check(found.size()==3,"the convergent integrals of a sector at one total power");
    // Massless lines: a squared line is an infrared divergence in four dimensions.
    auto momenta=family.momenta;
    for(auto& line:momenta.lines)line.mass_squared=Rational(0);
    const ibp::PropagatorBasis massless(ibp::quadratic_family(momenta,x,family.physical_count));
    const fm::Convergence soft(massless,1,Rational(4)),six(massless,1,Rational(6));
    check(!soft(integral(2,1)) && !soft(integral(1,2)) && !soft(integral(2,2)),"massless bubble in four dimensions");
    check(six(integral(2,2)) && !six(integral(1,3)),"massless bubble in six dimensions");
    // Exchange on a closure of this family in four dimensions: masters
    // [1,0] and [1,1] fail, the closure holds [3,0] and [2,1].
    level::Result level;level.success=true;
    level.ordered_basis={integral(1,0),integral(1,1)};
    level.matrix=Matrix{{e("eps/x"),e("0")},{e("1/(x-1)"),e("2*eps/(x-1)")}};
    level.requested_integrals={integral(2,1),integral(3,0)};
    level.target_rows=Matrix{{e("x"),e("eps+1")},{e("(eps+2)/x"),e("0")}};
    std::vector<fm::Step> steps;
    const auto exchanged=fm::exchange(level,0,1,four,steps,1);
    check(steps.size()==2 && exchanged.ordered_basis==std::vector<ibp::Integral>{integral(3,0),integral(2,1)},"both masters exchanged");
    check(exchanged.target_rows.size()==1 && exchanged.target_rows[0][0].is_zero() && exchanged.target_rows[0][1].is_one(),"the source is a master of the new basis");
    check(exchanged.requested_integrals.size()==1 && exchanged.source_identities.empty() && exchanged.differential_witnesses.empty(),"rows for the sources only");
    // The same steps from a record.
    std::vector<fm::Step> again;
    const auto replayed=fm::exchange(level,0,1,four,again,1,&steps);
    check(replayed.matrix==exchanged.matrix && replayed.ordered_basis==exchanged.ordered_basis,"recorded steps give the same closure");
    // A reduction with a moving pole the connection does not have is not
    // used; one with a new singular point in x alone is.
    auto pole=level;pole.target_rows[0]={e("x"),e("(eps+1)/(x-3+eps)")};
    fm::exchange(pole,0,1,four,steps,1);
    check(steps.size()==1 && steps[0].to==integral(3,0),"no step that adds a moving pole");
    const auto unguarded=fm::exchange(pole,0,1,four,steps,1,nullptr,false);
    check(steps.size()==2 && !fm::regular(pole.matrix,unguarded.matrix,0,1) && fm::regular(pole.matrix,exchanged.matrix,0,1),"the finished basis is checked instead");
    pole.target_rows[0]={e("x"),e("(eps+1)/(x-3)")};
    fm::exchange(pole,0,1,four,steps,1);
    check(steps.size()==2,"an apparent singular point in x is admitted");
    // Nor one whose connection has an epsilon pole no shear removes: J=[2,1]
    // defined with a coefficient that vanishes at epsilon = 0 at its master.
    auto dependent=level;dependent.target_rows[0]={e("1"),e("eps")};
    dependent.matrix=Matrix{{e("eps/x"),e("1/x")},{e("1/(x-1)"),e("2*eps/(x-1)")}};
    fm::exchange(dependent,0,1,four,steps,1);
    check(std::none_of(steps.begin(),steps.end(),[&](const fm::Step& step){return step.to==integral(2,1);}),"no step to masters that are dependent at epsilon zero");
    // choose(): candidates are asked for, all masters or none, and the
    // choice is recorded.
    const auto directory=std::filesystem::temp_directory_path()/("diffexp-finite-masters-"+std::to_string(::getpid()));
    std::filesystem::remove_all(directory);
    auto bare=level;bare.requested_integrals={integral(2,1)};bare.target_rows={level.target_rows[0]};
    std::size_t reductions=0;
    const auto reduce=[&](const std::vector<ibp::Integral>& requests,double) {
      ++reductions;
      auto wider=bare;wider.requested_integrals=requests;wider.target_rows.clear();
      for(const auto& request:requests) {
        if(request==integral(2,1))wider.target_rows.push_back(level.target_rows[0]);
        else if(request==integral(3,0))wider.target_rows.push_back(level.target_rows[1]);
        else if(request[1]==0)wider.target_rows.push_back({e("1/x^2"),e("0")});
        else wider.target_rows.push_back({e("x+eps"),e("x+2+eps")});
      }
      return wider;
    };
    fm::Report report;
    const auto chosen=fm::choose(massive,1,Rational(4),bare,{integral(2,1)},0,1,reduce,report,directory);
    check(report.divergent_before==2 && report.divergent_after==0 && report.steps.size()==2 && reductions==1 && !report.recorded,"choice with one further reduction");
    check(chosen.ordered_basis==std::vector<ibp::Integral>{integral(3,0),integral(2,1)} && chosen.target_rows.size()==1,"chosen masters");
    fm::Report second;
    const auto repeated=fm::choose(massive,1,Rational(4),bare,{integral(2,1)},0,1,reduce,second,directory);
    check(second.recorded && reductions==2 && repeated.matrix==chosen.matrix && repeated.ordered_basis==chosen.ordered_basis,"recorded choice repeated");
    // A divergent request keeps the masters, without a reduction.
    fm::Report kept;
    auto divergent=bare;divergent.requested_integrals={integral(1,1)};
    const auto same=fm::choose(massive,1,Rational(4),divergent,{integral(1,1)},0,1,reduce,kept,directory);
    check(kept.steps.empty() && reductions==2 && same.ordered_basis==bare.ordered_basis && !kept.note.empty(),"divergent requests keep the masters");
    // No finite replacement for one master: all masters are kept.
    const auto partial=[&](const std::vector<ibp::Integral>& requests,double seconds) {
      auto wider=reduce(requests,seconds);
      for(std::size_t i=0;i<requests.size();++i)if(requests[i][1]==0)wider.target_rows[i]={e("0"),e("0")};
      return wider;
    };
    fm::Report none;
    const auto unchanged=fm::choose(massive,1,Rational(4),bare,{integral(2,1)},0,1,partial,none);
    check(none.steps.empty() && none.divergent_after==2 && unchanged.ordered_basis==bare.ordered_basis && unchanged.matrix==bare.matrix,"all masters or none");
    std::filesystem::remove_all(directory);
  }
  std::cout<<"finite master tests passed\n";
}catch(const std::exception& failure){std::cerr<<failure.what()<<'\n';return 1;}}

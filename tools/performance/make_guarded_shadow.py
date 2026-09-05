#!/usr/bin/env python3
"""Generate an isolated adjoint-header candidate and a reviewable patch."""
from pathlib import Path
import difflib
root=Path(__file__).resolve().parents[2]
p=root/'include/diffexp/adjoint_transport.hpp';original=p.read_text()
if 'inline Boundary scalar_chart(' in original:
    start=original.index('inline Boundary scalar_chart(')
    split=original.index('struct ChartAcbArray',start)
    end=original.index('\ninline bool needs_centered_action',split)
    baseline=original[:start]+original[start:split].replace('inline Boundary scalar_chart(', 'inline Boundary chart(')+original[end:]
    for mode,content in [('baseline',baseline),('shadow',original)]:
        target=root/f'tools/performance/{mode}/diffexp/adjoint_transport.hpp'
        target.parent.mkdir(parents=True,exist_ok=True);target.write_text(content)
    patch=''.join(difflib.unified_diff(baseline.splitlines(True),original.splitlines(True),fromfile='a/include/diffexp/adjoint_transport.hpp',tofile='b/include/diffexp/adjoint_transport.hpp',n=0))
    (root/'tools/performance/guarded-dot.patch').write_text(patch)
    raise SystemExit(0)
start=original.index('inline Boundary chart(');end=original.index('\ninline bool needs_centered_action',start)
scalar=original[start:end].replace('inline Boundary chart(', 'inline Boundary scalar_chart(')
kernels=(root/'tools/performance/dot_kernels.hpp').read_text()
a=kernels.index('struct AcbArray');b=kernels.index('\nvoid require(',a)
dot=kernels[a:b].replace('AcbArray','ChartAcbArray').replace('Boundary experimental_chart(', 'inline Boundary dot_chart(')
dot=dot.replace('unsigned order,unsigned grouped)', 'unsigned order)')
a=dot.index('        if(grouped){');b=dot.index('\n      }',a)
dot=dot[:a]+'''        acb_dot(sum.raw(),target,0,coefficient,1,
            at(n,column,k-entries[e].epsilon),-stride,n+1,bits);
        acb_swap(target,sum.raw());'''+dot[b:]
wrapper='''
// Apply the same existing arithmetic reserve and per-chart radius guard.
// If grouped rounding consumes the reserve, retain the scalar reference path.
inline Boundary chart(const std::vector<Entry>& entries,const Boundary& initial,
    const Jet::Ball& center,const Jet::Ball& step,unsigned order) {
  auto output=dot_chart(entries,initial,center,step,order);
  if(needs_rational_cross_check(initial,output))
    return scalar_chart(entries,initial,center,step,order);
  return output;
}
'''
# Ball alias is deliberately local to the candidate namespace.
dot=dot.replace('const B& center,const B& step', 'const Jet::Ball& center,const Jet::Ball& step').replace('  const unsigned d=initial.size()', '  using B=Jet::Ball;\n  const unsigned d=initial.size()')
updated=original[:start]+scalar+'\n'+dot+wrapper+original[end:]
target=root/'tools/performance/shadow/diffexp/adjoint_transport.hpp';target.parent.mkdir(parents=True,exist_ok=True);target.write_text(updated)
patch=''.join(difflib.unified_diff(original.splitlines(True),updated.splitlines(True),fromfile='a/include/diffexp/adjoint_transport.hpp',tofile='b/include/diffexp/adjoint_transport.hpp',n=0))
(root/'tools/performance/guarded-dot.patch').write_text(patch)

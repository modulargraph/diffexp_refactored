// Included by rational_circuit.hpp inside namespace diffexp::rational_circuit.
//
// The Taylor recurrence of a chart on compact numbers. A chart performs
// millions of multiply-adds of an exact short coefficient (a polynomial
// coefficient shifted to a binary centre, mostly one limb) with a
// working-precision value. On ball structures each of them costs a cache miss
// or two, because every value is a 96-byte struct pointing at two heap
// blocks; here a value is 2L contiguous limbs, the orders a product sum walks
// through are adjacent in memory, and a product is one mpn_addmul_1.
//
// Number format. A value is (+-R + i*(+-I)) * 2^(64*(e-L)) with integers
// R,I < 2^(64L) and one exponent e for both parts, counted in limbs: no bit
// shift is needed to align a product. L is one limb more than the precision,
// for the up to 63 leading zero bits this alignment allows. One radius
// rho*2^(64*r), a double and a limb exponent, bounds the error of each part:
// the arithmetic is midpoint-radius ball arithmetic with the semantics of
// expand_balls, including its majorant bound for quotient streams.
namespace compact {
inline constexpr slong zero_exponent=std::numeric_limits<slong>::min()/4;
// m*2^(64e): radii and sizes. Operations round up by a relative 2^-40.
struct Scaled {
  double m=0;
  slong e=0;
};
inline double power64(slong k) {
  // 2^(64k); below the table the smallest entry still bounds from above.
  static const auto table=[]{std::array<double,31> t{};for(int i=0;i<31;++i)t[i]=std::ldexp(1.0,64*(i-15));return t;}();
  if(k>15)return std::numeric_limits<double>::infinity();
  return table[k< -15?0:k+15];
}
inline constexpr double up=1.0+0x1p-40;
inline Scaled normalized(Scaled x) {
  if(!(x.m>0))return {0,0};
  if(x.m<0x1p-200 || x.m>0x1p200) {
    if(!std::isfinite(x.m))return x;
    int binary;const double fraction=std::frexp(x.m,&binary);
    const slong whole=binary>=0?binary/64:-((-binary+63)/64);
    x.m=std::ldexp(fraction,static_cast<int>(binary-64*whole));x.e+=whole;
  }
  return x;
}
inline Scaled add(Scaled a,Scaled b) {
  if(!(a.m>0))return b;
  if(!(b.m>0))return a;
  if(a.e<b.e)std::swap(a,b);
  return normalized({(a.m+b.m*power64(b.e-a.e))*up,a.e});
}
inline Scaled multiply(Scaled a,Scaled b){return normalized({a.m*b.m*up,a.e+b.e});}
inline Scaled from_mag(const mag_t x) {
  if(mag_is_zero(x))return {0,0};
  if(mag_is_inf(x))return {std::numeric_limits<double>::infinity(),0};
  // x = man*2^(exp-30) with a 30-bit mantissa
  const slong binary=MAG_EXP(x);
  const slong whole=binary>=0?binary/64:-((-binary+63)/64);
  return {std::ldexp(static_cast<double>(MAG_MAN(x)),static_cast<int>(binary-64*whole-MAG_BITS))*up,whole};
}
inline void to_mag(mag_t out,Scaled x) {
  if(!(x.m>0)){mag_zero(out);return;}
  if(!std::isfinite(x.m)){mag_inf(out);return;}
  mag_set_d(out,x.m*up);mag_mul_2exp_si(out,out,64*x.e);
}
// +-limb*2^(64*exponent), an integer of `count` limbs (0: zero). One limb is
// held in place; a longer one lives in the Arena of its chart.
struct Short {
  static constexpr int capacity=9;
  const mp_limb_t* limb=&first;
  mp_limb_t first=0;
  int count=0;
  int exponent=0;
  bool negative=false;
  Short()=default;
  Short(const Short& other){*this=other;}
  Short& operator=(const Short& other) {
    first=other.first;count=other.count;exponent=other.exponent;negative=other.negative;
    limb=other.limb==&other.first?&first:other.limb;return *this;
  }
};
using Arena=std::deque<std::array<mp_limb_t,Short::capacity>>;
struct Coefficient {
  Short real,imaginary;
  slong top=zero_exponent;   // limbs: |part| < 2^(64*top)
  double size=0;             // (|real|+|imaginary|)/2^(64*top) or more
  double radius=0;           // radius of each part, in the same unit
};
inline bool shorten(Short& out,const arf_t x,Arena& arena) {
  out=Short{};
  if(arf_is_zero(x))return true;
  if(arf_is_special(x) || COEFF_IS_MPZ(*ARF_EXPREF(x)))return false;
  mp_srcptr digits;mp_size_t n;ARF_GET_MPN_READONLY(digits,n,x);
  if(n>Short::capacity-1)return false;
  mp_limb_t m[Short::capacity+1]={};for(mp_size_t i=0;i<n;++i)m[i]=digits[i];
  slong bit=ARF_EXP(x)-64*n;                 // value = m*2^bit
  while(n>1 && !m[0]){for(mp_size_t i=0;i+1<n;++i)m[i]=m[i+1];m[--n]=0;bit+=64;}
  const int low=static_cast<int>(flint_ctz(m[0]));
  if(low){mpn_rshift(m,m,n,low);bit+=low;if(n>1 && !m[n-1])--n;}
  const slong exponent=bit>=0?bit/64:-((-bit+63)/64);
  const int shift=static_cast<int>(bit-64*exponent);
  if(shift){m[n]=mpn_lshift(m,m,n,shift);if(m[n])++n;}
  if(n>Short::capacity || exponent>1000000000 || exponent< -1000000000)return false;
  out.count=static_cast<int>(n);out.first=m[0];
  if(n>1){arena.emplace_back();for(mp_size_t i=0;i<n;++i)arena.back()[i]=m[i];out.limb=arena.back().data();}
  out.exponent=static_cast<int>(exponent);out.negative=ARF_SGNBIT(x)!=0;
  return true;
}
inline bool shorten(Coefficient& out,const acb_t z,Arena& arena) {
  if(!acb_is_finite(z))return false;
  if(!shorten(out.real,arb_midref(acb_realref(z)),arena) || !shorten(out.imaginary,arb_midref(acb_imagref(z)),arena))return false;
  out.top=zero_exponent;out.size=0;out.radius=0;
  for(const auto* part:{&out.real,&out.imaginary})if(part->count)out.top=std::max<slong>(out.top,slong(part->exponent)+part->count);
  mag_t a,b;mag_init(a);mag_init(b);
  mag_max(a,arb_radref(acb_realref(z)),arb_radref(acb_imagref(z)));
  if(!mag_is_zero(a)) {
    // An inexact coefficient: its radius needs a unit even if the midpoint is zero.
    const auto r=from_mag(a);
    if(out.top==zero_exponent)out.top=r.e+1;
    out.radius=r.m*power64(r.e-out.top)*up;
    if(!std::isfinite(out.radius)){mag_clear(a);mag_clear(b);return false;}
  }
  if(out.top!=zero_exponent) {
    arf_get_mag(a,arb_midref(acb_realref(z)));arf_get_mag(b,arb_midref(acb_imagref(z)));mag_add(a,a,b);
    const auto s=from_mag(a);out.size=s.m*power64(s.e-out.top)*up;
  }
  mag_clear(a);mag_clear(b);
  return true;
}
// Values of one chart: 2L limbs each; exponents, flags and radii apart.
class Pool {
 public:
  enum : unsigned char { real_negative=1, imaginary_negative=2, real_zero=4, imaginary_zero=8 };
  Pool(slong limbs,std::size_t count)
      : limbs_(limbs),mantissa_(new mp_limb_t[count*2*static_cast<std::size_t>(limbs)]),exponent_(count,zero_exponent),
        flags_(count,real_zero|imaginary_zero),radius_(count) {}
  slong limbs()const{return limbs_;}
  mp_ptr real(std::size_t i){return mantissa_.get()+i*2*limbs_;}
  mp_ptr imaginary(std::size_t i){return real(i)+limbs_;}
  slong& exponent(std::size_t i){return exponent_[i];}
  unsigned char& flags(std::size_t i){return flags_[i];}
  Scaled& radius(std::size_t i){return radius_[i];}
  bool zero(std::size_t i)const{return exponent_[i]==zero_exponent;}
  void clear(std::size_t i){exponent_[i]=zero_exponent;flags_[i]=real_zero|imaginary_zero;radius_[i]={0,0};}
  void copy(std::size_t to,std::size_t from) {
    if(to==from)return;
    if(exponent_[from]!=zero_exponent)std::memcpy(real(to),real(from),sizeof(mp_limb_t)*2*limbs_);
    exponent_[to]=exponent_[from];flags_[to]=flags_[from];radius_[to]=radius_[from];
  }
  // Rounds toward zero; the lost part enters the radius.
  void store(std::size_t i,const acb_t z) {
    const arf_struct* parts[2]={arb_midref(acb_realref(z)),arb_midref(acb_imagref(z))};
    slong e=zero_exponent;
    for(const auto* part:parts)if(!arf_is_zero(part)) {
      const slong bits=ARF_EXP(part);
      e=std::max(e,bits>0?(bits+63)/64:-((-bits)/64));
    }
    mag_t r;mag_init(r);mag_max(r,arb_radref(acb_realref(z)),arb_radref(acb_imagref(z)));
    radius_[i]=from_mag(r);mag_clear(r);
    unsigned char flags=0;
    if(e==zero_exponent){exponent_[i]=e;flags_[i]=real_zero|imaginary_zero;return;}
    bool lost=false;
    for(int which=0;which<2;++which) {
      mp_ptr out=which?imaginary(i):real(i);
      for(slong k=0;k<limbs_;++k)out[k]=0;
      const arf_struct* part=parts[which];
      if(arf_is_zero(part)){flags|=which?imaginary_zero:real_zero;continue;}
      if(ARF_SGNBIT(part))flags|=which?imaginary_negative:real_negative;
      mp_srcptr digits;mp_size_t n;ARF_GET_MPN_READONLY(digits,n,part);
      // digits*2^(ARF_EXP-64n) = out*2^(64(e-L))
      const slong shift=ARF_EXP(part)-64*n-64*(e-limbs_);
      const slong whole=shift>=0?shift/64:-((-shift+63)/64);const int bits=static_cast<int>(shift-64*whole);
      mp_limb_t carry=0;bool any=false;
      for(mp_size_t k=0;k<=n;++k) {
        const mp_limb_t limb=k<n?digits[k]:0;
        const mp_limb_t value=bits?(limb<<bits)|carry:limb;
        if(bits)carry=limb>>(64-bits);
        const slong at=whole+k;
        if(at<0){if(value)lost=true;continue;}
        if(at<limbs_){out[at]=value;if(value)any=true;}
      }
      if(!any)flags|=which?imaginary_zero:real_zero;
    }
    exponent_[i]=e;flags_[i]=flags;
    if(lost)radius_[i]=add(radius_[i],{1.0,e-limbs_});
  }
  void load(acb_t z,std::size_t i) {
    mag_t r;mag_init(r);to_mag(r,radius_[i]);
    for(int which=0;which<2;++which) {
      arb_ptr out=which?acb_imagref(z):acb_realref(z);
      const bool none=zero(i) || (flags_[i]&(which?imaginary_zero:real_zero));
      mp_srcptr digits=which?imaginary(i):real(i);
      mp_size_t n=none?0:limbs_;while(n && !digits[n-1])--n;
      if(!n)arf_zero(arb_midref(out));
      else {
        arf_set_mpn(arb_midref(out),digits,n,(flags_[i]&(which?imaginary_negative:real_negative))!=0);
        arf_mul_2exp_si(arb_midref(out),arb_midref(out),64*(exponent_[i]-limbs_));
      }
      mag_set(arb_radref(out),r);
    }
    mag_clear(r);
  }
  // (|real|+|imaginary|)/2^(64e) or more, from the leading limbs.
  double fraction(std::size_t i) {
    const mp_limb_t top=std::max(real(i)[limbs_-1],imaginary(i)[limbs_-1]);
    return 2.0*(static_cast<double>(top)+2.0)*0x1p-64*up;
  }
  // log2 of the modulus, to the accuracy the tail estimate needs.
  double log2_size(std::size_t i) {
    if(zero(i))return -std::numeric_limits<double>::infinity();
    const auto leading=[&](mp_srcptr part) {
      return static_cast<double>(part[limbs_-1])+(limbs_>1?static_cast<double>(part[limbs_-2])*0x1p-64:0.0);
    };
    const double a=(flags_[i]&real_zero)?0.0:leading(real(i)),b=(flags_[i]&imaginary_zero)?0.0:leading(imaginary(i));
    const double size=std::hypot(a,b);
    if(!(size>0))return 64.0*(exponent_[i]-2); // only limbs below the two leading ones
    return std::log2(size)+64.0*(exponent_[i]-1);
  }
 private:
  slong limbs_;
  std::unique_ptr<mp_limb_t[]> mantissa_;
  std::vector<slong> exponent_;
  std::vector<unsigned char> flags_;
  std::vector<Scaled> radius_;
};
// value(top-a) is multiplied by coefficient[a], a<count. `midpoint` leaves
// the radii of the values out (the caller bounds them otherwise).
struct Line {
  const Coefficient* coefficient;
  unsigned count;
  std::size_t top;
  bool midpoint;
};
// out = the `limbs` limbs of sum from limb `low` up (zeros below limb zero);
// returns whether any is nonzero and collects the limbs dropped below.
inline bool leading(mp_ptr out,mp_srcptr sum,slong low,slong limbs,mp_limb_t& dropped) {
  mp_limb_t any=0;
  if(low>=0) {
    for(slong k=0;k<limbs;++k){out[k]=sum[low+k];any|=out[k];}
    for(slong at=0;at<low;++at)dropped|=sum[at];
  } else {
    for(slong k=0;k<-low;++k)out[k]=0;
    for(slong k=-low;k<limbs;++k){out[k]=sum[low+k];any|=out[k];}
  }
  return any!=0;
}
// Fixed-point sums of products below 2^(64*base): the integer X of `size`
// limbs, two's complement, stands for X*2^(64*(base+1-size)). The top limb
// takes carries and the sign; a product pass whose low limbs fall below limb
// zero loses less than one unit of limb one and is counted.
class Accumulator {
 public:
  static constexpr slong capacity=80;
  void start(slong size,slong base) {
    size_=size;base_=base;truncated_=0;
    for(slong k=0;k<size;++k)real_[k]=imaginary_[k]=0;
  }
  slong base()const{return base_;}
  slong size()const{return size_;}
  slong truncated()const{return truncated_;}
  // sum[at..) += or -= source*limb, carried into the limbs above.
  inline __attribute__((always_inline)) void add(mp_ptr sum,mp_limb_t limb,slong at,mp_srcptr source,slong length,bool negative) {
    if(__builtin_expect(at<0,0)) {
      ++truncated_;
      if(-at>=length)return;
      source-=at;length+=at;at=0;
    }
    mp_ptr high=sum+at+length;
    if(!negative) {
      const mp_limb_t carry=mpn_addmul_1(sum+at,source,length,limb);
      if(__builtin_expect((*high+=carry)<carry,0))for(mp_ptr q=high+1;q<sum+size_ && !++*q;++q){}
    } else {
      const mp_limb_t borrow=mpn_submul_1(sum+at,source,length,limb);
      const mp_limb_t before=*high;*high=before-borrow;
      if(__builtin_expect(before<borrow,0))for(mp_ptr q=high+1;q<sum+size_ && !(*q)--;++q){}
    }
  }
  inline __attribute__((always_inline)) void add(bool imaginary,const Short& c,mp_srcptr value,slong limbs,slong exponent,bool negative) {
    mp_ptr sum=imaginary?imaginary_:real_;
    const slong offset=slong(c.exponent)+exponent-limbs-base_-1+size_;
    if(__builtin_expect(c.count==1,1)){add(sum,c.first,offset,value,limbs,negative);return;}
    for(int j=0;j<c.count;++j)if(const mp_limb_t limb=c.limb[j])add(sum,limb,offset+j,value,limbs,negative);
  }
  // Writes the sum as value `target`; `radius` is the propagated radius in
  // units of 2^(64*base).
  void finish(Pool& pool,std::size_t target,double radius) {
    const slong limbs=pool.limbs();
    unsigned char flags=0;slong tops[2];
    for(int which=0;which<2;++which) {
      mp_ptr sum=which?imaginary_:real_;
      if(sum[size_-1]>>63){mpn_neg(sum,sum,size_);flags|=which?Pool::imaginary_negative:Pool::real_negative;}
      slong top=size_-1;while(top>=0 && !sum[top])--top;
      tops[which]=top;
    }
    const slong top=std::max(tops[0],tops[1]);
    Scaled total{radius,base_};
    if(truncated_)total=compact::add(total,Scaled{static_cast<double>(truncated_),base_+2-size_});
    if(top<0){pool.clear(target);pool.radius(target)=normalized(total);return;}
    const slong low=top-limbs+1,exponent=top+2+base_-size_;
    mp_limb_t dropped=0;
    for(int which=0;which<2;++which) {
      mp_srcptr sum=which?imaginary_:real_;mp_ptr out=which?pool.imaginary(target):pool.real(target);
      if(!leading(out,sum,low,limbs,dropped))flags|=which?Pool::imaginary_zero:Pool::real_zero;
    }
    pool.exponent(target)=exponent;pool.flags(target)=flags;
    // Exact sums stay exact: only limbs that held something cost a unit.
    if(dropped)total=compact::add(total,Scaled{1.0,exponent-limbs});
    pool.radius(target)=normalized(total);
  }
 private:
  mp_limb_t real_[capacity],imaginary_[capacity];
  slong size_=0,base_=0,truncated_=0;
};
// A full-precision factor (an inverse pivot).
struct Wide {
  std::vector<mp_limb_t> real,imaginary;
  slong exponent=zero_exponent;
  unsigned char flags=Pool::real_zero|Pool::imaginary_zero;
  Scaled radius,size;
};
struct Scratch {
  Accumulator sum;
  std::vector<Line> lines;
  std::vector<mp_limb_t> product,other,quotient;
  explicit Scratch(slong limbs):product(2*limbs+2),other(2*limbs+2),quotient(2*limbs+2) {}
};
inline constexpr std::size_t npos=std::numeric_limits<std::size_t>::max();
// target = initial +- sum over the lines, where `initial` is a value index or
// npos for zero. Returns the number of products. `exact_initial` leaves the
// radius of `initial` out, as Line::midpoint does for a line.
inline std::size_t dot(Pool& pool,Scratch& scratch,std::size_t target,std::size_t initial,bool subtract,
    bool exact_initial=false) {
  static const Short one=[]{Short s;s.first=1;s.count=1;return s;}();
  const slong limbs=pool.limbs();
  slong base=zero_exponent;
  if(initial!=npos && !pool.zero(initial))base=pool.exponent(initial)+1;
  for(const auto& line:scratch.lines)for(unsigned a=0;a<line.count;++a) {
    const slong top=line.coefficient[a].top,e=pool.exponent(line.top-a);
    if(top!=zero_exponent && e!=zero_exponent && top+e>base)base=top+e;
  }
  // Radii, in units of 2^(64*base); with no midpoint at all, of the first radius.
  slong unit=base;double radius=0;
  const auto widen=[&](double m,slong e) {
    if(!(m>0))return;
    if(unit==zero_exponent)unit=e;
    radius+=m*power64(e-unit);
  };
  if(initial!=npos && !exact_initial)widen(pool.radius(initial).m,pool.radius(initial).e);
  std::size_t terms=0;
  auto& sum=scratch.sum;
  if(base!=zero_exponent) {
    sum.start(limbs+5,base);
    if(initial!=npos && !pool.zero(initial)) {
      const auto flags=pool.flags(initial);const slong e=pool.exponent(initial);
      if(!(flags&Pool::real_zero))sum.add(false,one,pool.real(initial),limbs,e,(flags&Pool::real_negative)!=0);
      if(!(flags&Pool::imaginary_zero))sum.add(true,one,pool.imaginary(initial),limbs,e,(flags&Pool::imaginary_negative)!=0);
    }
  }
  for(const auto& line:scratch.lines)for(unsigned a=0;a<line.count;++a) {
    const auto& c=line.coefficient[a];const std::size_t index=line.top-a;
    if(c.top==zero_exponent)continue;
    const auto& r=pool.radius(index);
    if(!line.midpoint && r.m>0)widen(c.size*r.m,c.top+r.e);
    if(pool.zero(index))continue;
    const auto flags=pool.flags(index);const slong e=pool.exponent(index);
    if(c.radius>0)widen(c.radius*(pool.fraction(index)+r.m*power64(r.e-e)),c.top+e);
    ++terms;
    const bool vr=!(flags&Pool::real_zero),vi=!(flags&Pool::imaginary_zero);
    const bool nr=(flags&Pool::real_negative)!=0,ni=(flags&Pool::imaginary_negative)!=0;
    if(c.real.count) {
      if(vr)sum.add(false,c.real,pool.real(index),limbs,e,c.real.negative!=nr!=subtract);
      if(vi)sum.add(true,c.real,pool.imaginary(index),limbs,e,c.real.negative!=ni!=subtract);
    }
    if(c.imaginary.count) {
      if(vi)sum.add(false,c.imaginary,pool.imaginary(index),limbs,e,c.imaginary.negative!=ni!=!subtract);
      if(vr)sum.add(true,c.imaginary,pool.real(index),limbs,e,c.imaginary.negative!=nr!=subtract);
    }
  }
  if(base==zero_exponent) {
    pool.clear(target);pool.radius(target)=normalized({radius*up,unit==zero_exponent?0:unit});return 0;
  }
  sum.finish(pool,target,radius*up);
  return terms;
}
// value *= factor.
inline void multiply(Pool& pool,Scratch& scratch,std::size_t index,const Wide& factor) {
  const slong limbs=pool.limbs();
  // radius: |f| r + r_f |v| + r_f r, then the rounding of the product.
  const Scaled before=pool.radius(index);
  Scaled bound=multiply(factor.size,before);
  if(!pool.zero(index))bound=add(bound,multiply(factor.radius,{pool.fraction(index),pool.exponent(index)}));
  bound=add(bound,multiply(factor.radius,before));
  if(pool.zero(index) || factor.exponent==zero_exponent){pool.clear(index);pool.radius(index)=bound;return;}
  const auto flags=pool.flags(index);
  const slong size=2*limbs+2;
  // Both parts in two's complement over 2L+2 limbs.
  mp_limb_t* parts[2]={scratch.product.data(),scratch.other.data()};
  mp_ptr temporary=scratch.quotient.data();
  const auto accumulate=[&](mp_ptr sum,mp_srcptr a,mp_srcptr b,bool negative) {
    mpn_mul_n(temporary,a,b,limbs);
    if(!negative)mpn_add_1(sum+2*limbs,sum+2*limbs,2,mpn_add_n(sum,sum,temporary,2*limbs));
    else mpn_sub_1(sum+2*limbs,sum+2*limbs,2,mpn_sub_n(sum,sum,temporary,2*limbs));
  };
  for(auto* part:parts)for(slong k=0;k<size;++k)part[k]=0;
  const bool vr=!(flags&Pool::real_zero),vi=!(flags&Pool::imaginary_zero);
  const bool fr=!(factor.flags&Pool::real_zero),fi=!(factor.flags&Pool::imaginary_zero);
  const bool nvr=(flags&Pool::real_negative)!=0,nvi=(flags&Pool::imaginary_negative)!=0;
  const bool nfr=(factor.flags&Pool::real_negative)!=0,nfi=(factor.flags&Pool::imaginary_negative)!=0;
  if(vr&&fr)accumulate(parts[0],pool.real(index),factor.real.data(),nvr!=nfr);
  if(vi&&fi)accumulate(parts[0],pool.imaginary(index),factor.imaginary.data(),nvi!=nfi!=true);
  if(vr&&fi)accumulate(parts[1],pool.real(index),factor.imaginary.data(),nvr!=nfi);
  if(vi&&fr)accumulate(parts[1],pool.imaginary(index),factor.real.data(),nvi!=nfr);
  unsigned char out=0;slong tops[2];
  for(int which=0;which<2;++which) {
    mp_ptr sum=parts[which];
    if(sum[size-1]>>63){mpn_neg(sum,sum,size);out|=which?Pool::imaginary_negative:Pool::real_negative;}
    slong top=size-1;while(top>=0 && !sum[top])--top;
    tops[which]=top;
  }
  const slong top=std::max(tops[0],tops[1]);
  if(top<0){pool.clear(index);pool.radius(index)=bound;return;}
  const slong low=top-limbs+1,exponent=pool.exponent(index)+factor.exponent-2*limbs+top+1;
  mp_limb_t dropped=0;
  for(int which=0;which<2;++which) {
    if(!leading(which?pool.imaginary(index):pool.real(index),parts[which],low,limbs,dropped))out|=which?Pool::imaginary_zero:Pool::real_zero;
  }
  pool.exponent(index)=exponent;pool.flags(index)=out;
  if(dropped)bound=add(bound,{1.0,exponent-limbs});
  pool.radius(index)=bound;
}
// value /= divisor, a small positive integer.
inline void divide(Pool& pool,Scratch& scratch,std::size_t index,ulong divisor) {
  auto& radius=pool.radius(index);
  radius.m=radius.m/static_cast<double>(divisor)*up;
  if(pool.zero(index) || divisor==1)return;
  const slong limbs=pool.limbs();
  mp_ptr q[2]={scratch.quotient.data(),scratch.quotient.data()+limbs+1};
  auto flags=pool.flags(index);mp_limb_t high=0,inexact=0;
  for(int which=0;which<2;++which) {
    if(flags&(which?Pool::imaginary_zero:Pool::real_zero)){for(slong k=0;k<=limbs;++k)q[which][k]=0;continue;}
    inexact|=mpn_divrem_1(q[which],1,which?pool.imaginary(index):pool.real(index),limbs,divisor);
    high|=q[which][limbs];
  }
  for(int which=0;which<2;++which) {
    mp_ptr to=which?pool.imaginary(index):pool.real(index);mp_limb_t any=0;
    for(slong k=0;k<limbs;++k){to[k]=q[which][k+(high?1:0)];any|=to[k];}
    if(high)inexact|=q[which][0];
    if(!any)flags|=which?Pool::imaginary_zero:Pool::real_zero;
  }
  if(!high)--pool.exponent(index);
  pool.flags(index)=flags;
  if(inexact)radius=add(radius,{1.0,pool.exponent(index)-limbs});
}
inline bool widen(Wide& out,const acb_t z,slong limbs) {
  for(const auto* part:{arb_midref(acb_realref(z)),arb_midref(acb_imagref(z))})
    if(!arf_is_zero(part) && (arf_is_special(part) || COEFF_IS_MPZ(*ARF_EXPREF(part))))return false;
  Pool single(limbs,1);single.store(0,z);
  out.real.assign(single.real(0),single.real(0)+limbs);out.imaginary.assign(single.imaginary(0),single.imaginary(0)+limbs);
  out.exponent=single.exponent(0);out.flags=single.flags(0);out.radius=single.radius(0);
  mag_t a,b;mag_init(a);mag_init(b);
  arb_get_mag(a,acb_realref(z));arb_get_mag(b,acb_imagref(z));mag_add(a,a,b);out.size=from_mag(a);
  mag_clear(a);mag_clear(b);
  return true;
}
// The series of a chart in compact form: y[(r*width+k)*stride+n], then
// scratch for the Horner values of evaluate().
struct Stored {
  Pool pool;
  std::size_t stride,horner;
  Stored(slong limbs,std::size_t count,std::size_t stride_,std::size_t horner_):pool(limbs,count),stride(stride_),horner(horner_) {}
};
} // namespace compact

// The recurrence of expand_balls on compact numbers. Returns nothing when a
// coefficient or value does not fit the format; the caller then uses
// expand_balls. `defect_radii` as there.
inline std::optional<Series> expand_compact(const Prepared& prepared,const Boundary& boundary,Statistics* statistics,
    bool defect_radii=false) {
  using namespace compact;
  const auto& p=prepared.data();const auto& plan=p.compiled.data();
  const unsigned width=unsigned(boundary[0].size()),d=plan.dimension,order=p.order;
  const slong limbs=(p.precision+63)/64+1;
  if(limbs+5>Accumulator::capacity)return std::nullopt;
  for(const auto& row:boundary)for(const auto& value:row)for(const auto* part:{acb_realref(value.raw()),acb_imagref(value.raw())})
    if(!arf_is_zero(arb_midref(part)) && COEFF_IS_MPZ(*ARF_EXPREF(arb_midref(part))))return std::nullopt;
  // Short coefficients and wide pivots of this centre.
  Arena arena;
  std::vector<std::vector<Coefficient>> polynomials(p.polynomials.size());
  for(std::size_t i=0;i<p.polynomials.size();++i) {
    polynomials[i].resize(p.polynomials[i].size());
    for(std::size_t a=0;a<p.polynomials[i].size();++a)
      if(!shorten(polynomials[i][a],p.polynomials[i][a].raw(),arena))return std::nullopt;
  }
  std::vector<Wide> pivots(p.inverse_pivots.size());
  for(std::size_t i=0;i<pivots.size();++i)if(!widen(pivots[i],p.inverse_pivots[i].raw(),limbs))return std::nullopt;
  // Structural reachability, as in expand_balls.
  std::vector<std::vector<unsigned>> dependents(d);
  for(const auto& entry:plan.entries)
    if(!entry.numerator.empty() && entry.numerator.front().epsilon<width)
      dependents[plan.streams[entry.stream].column].push_back(entry.row);
  std::vector<bool> active(d,false),active_stream(plan.streams.size(),false);
  Series series{d,width,order,p.precision,{},{},{}};
  auto& queue=series.rows;
  for(unsigned i=0;i<d;++i)for(const auto& value:boundary[i])if(!value.is_zero()){active[i]=true;queue.push_back(i);break;}
  for(std::size_t next=0;next<queue.size();++next)for(auto row:dependents[queue[next]])if(!active[row]){active[row]=true;queue.push_back(row);}
  std::vector<unsigned> row_index(d);
  for(unsigned i=0;i<queue.size();++i)row_index[queue[i]]=i;
  for(const auto& entry:plan.entries)
    if(active[plan.streams[entry.stream].column] && !entry.numerator.empty() && entry.numerator.front().epsilon<width)
      active_stream[entry.stream]=true;
  // Storage: y[(r*width+k)*(order+1)+n]; a quotient stream keeps a ring of R
  // orders per epsilon index, written twice (at s and s+R) so that the orders
  // n, n-1, ... are adjacent below position (n mod R)+R.
  const std::size_t stride=order+1,rows=queue.size();
  std::vector<std::size_t> offsets(plan.streams.size(),0);
  const std::size_t horner=rows*width*stride;
  std::size_t count=horner+rows*width;
  Statistics stats;stats.coefficient_cells=p.coefficient_cells;
  for(unsigned id=0;id<plan.streams.size();++id) {
    offsets[id]=count;
    if(active_stream[id] && p.ring_lengths[id]){count+=std::size_t(2)*p.ring_lengths[id]*width;stats.auxiliary_cells+=std::size_t(p.ring_lengths[id])*width;}
  }
  // The logical cell count of expand_balls, which budgets are stated in.
  stats.live_cells=std::size_t(order+1)*rows*width+stats.auxiliary_cells+std::size_t(d)*width+p.coefficient_cells+
    p.inverse_pivots.size()+2+2*p.dot_capacity+1;
  auto stored=std::make_shared<Stored>(limbs,count,stride,horner);
  Pool& pool=stored->pool;
  const auto yat=[&](unsigned n,unsigned row,unsigned k){return (std::size_t(row_index[row])*width+k)*stride+n;};
  // Index of Taylor order n of a stream; order n-a is at that index minus a.
  const auto line_top=[&](unsigned id,unsigned n,unsigned k)->std::size_t {
    const unsigned ring=p.ring_lengths[id];
    if(!ring)return yat(n,plan.streams[id].column,k);
    return offsets[id]+std::size_t(k)*2*ring+(n%ring)+ring;
  };
  for(auto i:queue)for(unsigned k=0;k<width;++k)pool.store(yat(0,i,k),boundary[i][k].raw());
  // Majorant state per active stream, as in expand_balls.
  std::vector<std::vector<Scaled>> chains(plan.streams.size());
  std::vector<std::vector<Scaled>> inverse_distances(p.inverse_distances.size());
  for(std::size_t i=0;i<p.inverse_distances.size();++i)
    for(std::size_t level=0;level<p.inverse_distances[i].size();++level)inverse_distances[i].push_back(from_mag(p.inverse_distances[i][level]));
  for(unsigned id=0;id<plan.streams.size();++id)
    if(active_stream[id] && p.ring_lengths[id] && plan.denominators[plan.streams[id].denominator].majorized)
      chains[id].assign(inverse_distances[plan.streams[id].denominator].size()*width,Scaled{});
  const unsigned stream_count=unsigned(plan.streams.size()),row_count=unsigned(rows);
  std::size_t order_work=0;
  for(unsigned id=0;id<stream_count;++id)if(active_stream[id]&&p.ring_lengths[id])
    for(auto b:plan.denominators[plan.streams[id].denominator].blocks)if(b.epsilon<width)
      order_work+=p.polynomials[b.polynomial].size()*(width-b.epsilon);
  for(const auto& entry:plan.entries)if(active_stream[entry.stream]&&active[entry.row])
    for(auto b:entry.numerator)if(b.epsilon<width)
      order_work+=p.polynomials[b.polynomial].size()*(width-b.epsilon);
  const bool threaded=order_work>=4*parallel_work;
  std::size_t dot_products=0,addmul_operations=0,scalar_operations=0;
#pragma omp parallel if(threaded) reduction(+:dot_products,addmul_operations,scalar_operations)
  {
    Scratch scratch(limbs);
    for(unsigned n=0;n<order;++n) {
#pragma omp for schedule(dynamic)
      for(unsigned id=0;id<stream_count;++id) {
        if(!active_stream[id] || !p.ring_lengths[id])continue;
        const auto& stream=plan.streams[id];
        const auto& denominator=plan.denominators[stream.denominator];
        const auto& inverse=inverse_distances[stream.denominator];
        const unsigned ring=p.ring_lengths[id];
        for(unsigned k=0;k<width;++k) {
          const std::size_t target=line_top(id,n,k);
          scratch.lines.clear();
          for(auto b:denominator.blocks)if(b.epsilon<=k) {
            const auto& q=polynomials[b.polynomial];
            const unsigned first=b.epsilon?0:1,last=std::min<unsigned>(unsigned(q.size())-1,n);
            if(first>last)continue;
            // Order n of a lower epsilon index of this stream was written
            // earlier in this loop; the orders below n are the ring.
            scratch.lines.push_back({q.data()+first,last-first+1,line_top(id,n,k-b.epsilon)-first,
              denominator.majorized && !b.epsilon});
          }
          addmul_operations+=dot(pool,scratch,target,yat(n,stream.column,k),true,defect_radii);++dot_products;
          multiply(pool,scratch,target,pivots[stream.denominator]);++scalar_operations;
          if(denominator.majorized) {
            // See expand_balls: the radius of the value bounds the defect of
            // the midpoint recurrence; the chain bounds the quotient.
            Scaled bound=pool.radius(target);bound.m*=2; // both parts
            for(std::size_t level=0;level<inverse.size();++level) {
              auto& previous=chains[id][level*width+k];
              bound=add(bound,multiply(previous,inverse[level]));
              previous=bound;
            }
            pool.radius(target)=bound;
          }
          pool.copy(target-ring,target);
        }
      }
#pragma omp for schedule(dynamic)
      for(unsigned r=0;r<row_count;++r)for(unsigned k=0;k<width;++k) {
        const auto row=queue[r];
        scratch.lines.clear();
        for(auto id:plan.row_entries[row]) {
          const auto& entry=plan.entries[id];
          if(!active_stream[entry.stream])continue;
          for(auto b:entry.numerator)if(b.epsilon<=k) {
            const auto& numerator=polynomials[b.polynomial];
            scratch.lines.push_back({numerator.data(),std::min<unsigned>(unsigned(numerator.size()),n+1),
              line_top(entry.stream,n,k-b.epsilon),defect_radii && !p.ring_lengths[entry.stream]});
          }
        }
        const std::size_t target=yat(n+1,row,k);
        addmul_operations+=dot(pool,scratch,target,npos,false);++dot_products;
        divide(pool,scratch,target,n+1);++scalar_operations;
      }
    }
  }
  stats.dot_products=dot_products;stats.addmul_operations=addmul_operations;stats.scalar_operations=scalar_operations;
  series.compact=std::move(stored);
  if(statistics)*statistics=stats;
  return series;
}
// Ball coefficients of a compact series, for consumers that read them.
inline void materialize(Series& series) {
  if(!series.compact || !series.y.empty())return;
  auto& stored=*std::static_pointer_cast<compact::Stored>(series.compact);
  const std::size_t rows=series.rows.size();
  series.y.resize(std::size_t(series.order+1)*rows*series.width);
  for(unsigned n=0;n<=series.order;++n)for(std::size_t r=0;r<rows;++r)for(unsigned k=0;k<series.width;++k)
    stored.pool.load(series.y[(std::size_t(n)*rows+r)*series.width+k].raw(),(r*series.width+k)*stored.stride+n);
}
// Horner value and tail estimate of a compact series; nothing when the step
// is not a short exact number.
inline std::optional<Boundary> evaluate_compact(const Series& series,const B& step,double* truncation_log2) {
  using namespace compact;
  Coefficient factor;Arena arena;
  if(!shorten(factor,step.raw(),arena) || factor.radius>0 || factor.top==zero_exponent)return std::nullopt;
  auto& stored=*std::static_pointer_cast<Stored>(series.compact);
  Pool& pool=stored.pool;
  Scratch scratch(pool.limbs());
  Boundary result;result.reserve(series.dimension);
  for(unsigned i=0;i<series.dimension;++i)result.emplace_back(series.width);
  const double step_log2=TaylorTailEstimate::log2_midpoint(step.raw());
  const unsigned first=series.order>3?series.order-3:1;
  double worst=-std::numeric_limits<double>::infinity();
  for(std::size_t r=0;r<series.rows.size();++r)for(unsigned k=0;k<series.width;++k) {
    const std::size_t base=(r*series.width+k)*stored.stride,value=stored.horner+r*series.width+k;
    pool.copy(value,base+series.order);
    for(unsigned n=series.order;n-->0;) {
      scratch.lines.assign(1,Line{&factor,1,value,false});
      dot(pool,scratch,value,base+n,false);
    }
    pool.load(result[series.rows[r]][k].raw(),value);
    if(!truncation_log2)continue;
    const double scale=pool.log2_size(value);
    double terms[4],top=-std::numeric_limits<double>::infinity();unsigned count=0;
    for(unsigned n=first;n<=series.order;++n) {
      const double term=pool.log2_size(base+n)+static_cast<double>(n)*step_log2;
      if(!std::isfinite(term))continue;
      terms[count++]=term;top=std::max(top,term);
    }
    if(!count)continue;
    double sum=0;for(unsigned i=0;i<count;++i)sum+=std::exp2(terms[i]-top);
    const double denominator=scale>0?scale+std::log2(1+std::exp2(-scale)):std::log2(1+std::exp2(scale));
    worst=std::max(worst,top+std::log2(sum)-denominator);
  }
  if(truncation_log2)*truncation_log2=worst;
  return result;
}

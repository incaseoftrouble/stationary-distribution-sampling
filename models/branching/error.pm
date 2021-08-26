dtmc

const double e;
const int m;
const double p;
module loop
  pos : [1 .. 2] init 1;
  x : [1 .. m] init 1;

  [] pos = 1 & x < m -> (1 - e) : (pos'=2) + e : (x'=min(x+1,m));
  [] pos = 1 & x = m -> (1 - e) : (pos'=2) + e : (x'=1);
  [] pos = 2 -> (1 - p - e) : (pos'=2) + p : (pos'=1) + e : (pos'=1) & (x'=max(x-1,1));
endmodule
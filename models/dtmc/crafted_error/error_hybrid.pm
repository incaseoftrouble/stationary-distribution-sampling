dtmc

const double e;
module loop
  pos : [1 .. 2] init 1;
  x : [1 .. 2] init 1;

  [] pos = 1 & x < 2 -> (1 - e) : (pos'=2) + e : (x'=2);
  [] pos = 1 & x = 2 -> (1 - e) : (pos'=2) + e : (x'=1);
  [] pos = 2         -> (1 - e) : (pos'=1) + e : (pos'=1) & (x'=1);
endmodule
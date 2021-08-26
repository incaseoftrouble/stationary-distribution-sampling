dtmc

const int size;

module x
  x : [1 .. size] init 1;
  [step] true -> 1/2 : (x'=min(x+1,size)) + (1 - 1/2) : (x'=1);
endmodule

module y
  y : [1 .. size] init 1;
  [step] true -> 1/2 : (y'=min(y+1,size)) + (1 - 1/2) : (y'=1);
endmodule

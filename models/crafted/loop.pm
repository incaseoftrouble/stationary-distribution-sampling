dtmc

const int size;

module x
  x : [1 .. size] init 1;
  [step] true -> (1 / 2 + x / (2 * size)) : (x'=max(x-1,1)) + (1 / 2 - x / (2 * size)) : (x'=min(x+1,size));
endmodule
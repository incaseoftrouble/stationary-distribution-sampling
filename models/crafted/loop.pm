dtmc

const int size;

module x
  x : [1 .. size] init 1;
  [step] true -> (x / size) : (x'=max(x-1,1)) + (1 - x / size) : (x'=min(x+1,size));
endmodule
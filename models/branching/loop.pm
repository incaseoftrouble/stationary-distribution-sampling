dtmc

const int size;

module loop
  loop : [1 .. size] init 1;

  [] loop < size -> (1/2) : (loop'=loop+1) + (1/2) : (loop'=loop);
  [] loop = size -> 1 : (loop'=1);
endmodule
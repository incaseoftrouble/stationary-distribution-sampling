dtmc

const int N;
const double p;

module state
  s : [1 .. N];

  [] true -> p : (s' = min(s+1,N)) + (1-p) : (s' = max(s-1,1));
endmodule


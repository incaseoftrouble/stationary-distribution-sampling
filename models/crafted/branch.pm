dtmc

const double p;
const int steps;
const int levels;
const int loops;

module level
  level : [1 .. levels] init 1;

  [step] step < steps -> p * 2 / 3 : (level'=min(level+1,levels))
                   + p * 1 / 3 : (level'=min(level*2,levels))
                   + (1-p) : (level'=level);
  [loop] step = steps -> 1 : true;
endmodule

module step
  step : [1 .. steps] init 1;

  [step] step < steps -> 1 : (step' = step+1);
  [loop] step = steps -> 1 : true;
endmodule

module loop
  loop : [1 .. loops] init 1;

  [step] true -> 1 : true;
  [loop] loop < loops -> (1/2) : (loop'=loop+1) + (1/2) : (loop'=1);
  [loop] loop = loops -> 1 : (loop'=1);
endmodule

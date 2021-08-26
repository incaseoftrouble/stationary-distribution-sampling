dtmc

const int steps;
const double p;
const int levels;

module level
  level : [1 .. levels] init 1;

  [step] step < steps -> p / 2 : (level'=min(level+1,levels))
                   + p / 4 : (level'=min(level+2,levels))
                   + p / 8 : (level'=min(level+3,levels))
                   + p / 8 : (level'=min(level*2,levels))
                   + (1-p) : (level'=level);
  [loop] step = steps -> 1 : true;
endmodule

module step
  step : [1 .. steps] init 1;

  [step] step < steps -> 1 : (step' = step+1);
  [loop] step = steps -> 1 : true;
endmodule
package stationary.util;

import de.tum.in.probmodels.util.Util;
import java.util.Locale;

public record Bound(double lower, double upper) {
  public static final Bound UNKNOWN = new Bound(0.0, 1.0);
  public static final Bound ONE = new Bound(1.0, 1.0);
  public static final Bound ZERO = new Bound(0.0, 0.0);

  public double difference() {
    return upper - lower;
  }

  public double average() {
    return (upper + lower) / 2.0;
  }

  @Override
  public String toString() {
    return String.format(Locale.ENGLISH, "[%.5f+%.5f]", lower, difference());
  }

  public boolean contains(Bound bound) {
    return Util.lessOrEqual(lower(), bound.lower()) && Util.lessOrEqual(bound.upper(), upper());
  }

  public static Bound of(double lower, double upper) {
    assert Util.lessOrEqual(lower, upper);
    return new Bound(lower, upper);
  }

  public static Bound exact(double value) {
    return new Bound(value, value);
  }
}

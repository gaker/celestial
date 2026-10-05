# Accuracy

The series shipped in the crate are truncated. They keep 83,780 of VSOP2013's
2,607,947 terms and 10,754 of ELP/MPP02's 35,901. Each body is cut where the
dropped terms move its direction, seen from the Earth, by at most 10 mas over
1950–2050, or by less than the full theory already misses JPL's ephemeris
where that is worse. 10 mas is about the aberration error that heliocentric
velocities carry anyway (see [Conventions](./theories.md#conventions)).

## Measured Against DE432s

Each body was compared with JPL's DE432s at 20,000 random TDB epochs from 1950
to 2050, the span DE432s covers. Both sides are geometric positions in ICRS.
The Earth comes from `Vsop2013Earth`.

As seen from the Earth's center:

| Body    | Direction, max (mas) | Direction, 99% (mas) | Distance, max (km) |
|---------|---------------------:|---------------------:|-------------------:|
| Sun     |                  2.3 |                  1.7 |                0.8 |
| Moon    |                   38 |                   30 |              0.067 |
| Mercury |                  7.4 |                  3.9 |                2.9 |
| Venus   |                  6.9 |                  4.1 |                1.3 |
| Mars    |                   15 |                  7.1 |                3.1 |
| Jupiter |                   30 |                   27 |                 33 |
| Saturn  |                  6.8 |                  5.0 |                 24 |
| Uranus  |                  752 |                  679 |               2942 |
| Neptune |                  176 |                  155 |               1542 |
| Pluto   |                 3041 |                 2804 |              70976 |

The Moon row is `ElpMpp02Moon::new()`, and its distance column is the largest
miss in position. `with_de405_fit()` gives 21 mas, 13 mas and 0.036 km.

The Sun's error is the Earth's: `Vsop2013Earth` is within 1.7 km of DE432s. As
a rough guide, Mercury to Saturn are good to 30 mas, Uranus and Neptune to
0.8″, and Pluto to 3″.

## Where the Error Comes From

The next table compares heliocentric errors. It sets the shipped tables beside
the full VSOP2013 series, measured against DE432s at 2,000 epochs over the
same span:

| Body                  | Shipped (km) | Shipped (mas) | Full series (mas) |
|-----------------------|-------------:|--------------:|------------------:|
| Mercury               |          3.3 |            14 |               4.9 |
| Venus                 |          1.0 |           1.8 |               1.0 |
| Earth-Moon barycenter |          1.7 |           2.3 |               1.0 |
| Mars                  |          4.5 |           4.3 |               3.2 |
| Jupiter               |          104 |            29 |                22 |
| Saturn                |           45 |           6.9 |               1.8 |
| Uranus                |        10354 |           749 |               562 |
| Neptune               |         4002 |           182 |               117 |
| Pluto                 |       115429 |          3122 |              2579 |

Mercury, Venus, Mars and Saturn are cut for the 10 mas limit. Mercury looks
worse here than from the Earth because it is so close to the Sun. Jupiter and
the planets beyond it are cut where the dropped terms are smaller than the
theory's own error, so keeping more terms would barely change their rows. The
Earth-Moon barycenter is cut more finely than either rule asks, because every
position seen from the Earth carries the Earth's error. Keeping more terms
costs memory and speed; see [Regenerating the Tables](./generators.md).

For the Moon, the full series misses DE432s by 60 m with `new()` and 21 m with
`with_de405_fit()`. The shipped table drops terms worth at most about 19 m.

## Outside 1950–2050

DE432s spans only 1950–2050, so the numbers above say nothing about other
dates. The dropped terms include secular ones (multiplied by powers of time),
so expect the error to grow with distance from J2000. The date ranges each
theory accepts (see [Planets, Sun and Moon](./theories.md#date-ranges)) are
far wider than the span measured here.

For anything tighter, or for barycentric positions, read a JPL kernel; see
[JPL SPK Kernels](./spk.md).

# Despirala-Engine

Finds the optimal strategy for the dice game Despirala.
Can also evaluate luck and mistakes made during play.

## Rules of the game:

1. At the start of each turn you roll all six dice for free and then gain five goods.
2. After that you choose to either reroll the dice (for the cost of one good) or choose some combination that you have not attempted before.
3. After picking a combination you have to complete it.
4. To that end, you repeatedly take some of the dice (of your choosing) and reroll only them (for the cost of one good).
5. If you successfully complete the combination, you win the amount of points it is worth and your turn ends.
6. If you run out of goods before completing the chosen combination, your turn ends and you do not get any points for that combination.
7. The six single number combinations (ones, twos, ...) are an exception: you choose when you are done rolling and get points proportional to the number of correct dice.
8. At the end of the game you get one bonus point per good you have left over.

## List of combinations:

| Name | Contents | Points |
| ---- | ---- | ---- |
| Collect ones | N ones | Sum of dice in the combination |
| Collect twos | N twos | Sum of dice in the combination |
| Collect threes | N threes | Sum of dice in the combination |
| Collect fours | N fours | Sum of dice in the combination |
| Collect fives | N fives | Sum of dice in the combination |
| Collect sixes | N sixes | Sum of dice in the combination |
| Three pairs | X, X, Y, Y, Z, Z | Sum of dice in the combination |
| Two triples | X, X, X, Y, Y, Y | Sum of dice in the combination |
| Four of a kind | X, X, X, X | 40 |
| Kamerun | 4, 5, 5, 6, 6, 6 | 45 |
| Straight | 1, 2, 3, 4, 5, 6 | 50 |
| Six of a kind | X, X, X, X, X, X | 60 |
| General | 6, 6, 6, 6, 6, 6 | 70 |
| Despirala | 1, 1, 1, 1, 1, 6 | 80 |

Note: for combinations such as pairs/triples/X of a kind, you have to pick which pairs/triples/X you are doing.
E.g. Three pairs: threes, fives and sixes.

## Performance stats

Stats for model_max are obtained with 40 million tests.
All othere stats are obtained with a million tests.
Note that for a perfect model, "expected score" will be equal to the empirical mean.

### Stats for model_max:

Expected score: 443.616

Mean: 443.617 \
Stdev: 61.4432 \
5th percentile: 310 \
25th percentile: 420 \
50th percentile: 468 \
75th percentile: 484 \
95th percentile: 501 \
Mode: 478

### Stats for model_1000:

Expected score: 443.694

Mean: 443.44 \
Stdev: 61.5872 \
5th percentile: 310 \
25th percentile: 420 \
50th percentile: 468 \
75th percentile: 483 \
95th percentile: 501 \
Mode: 477

### Stats for model_50:

Expected score: 444.524

Mean: 442.368 \
Stdev: 62.1897 \
5th percentile: 308 \
25th percentile: 418 \
50th percentile: 467 \
75th percentile: 483 \
95th percentile: 500 \
Mode: 480

### Stats for model_5:

Expected score: 457.534

Mean: 425.704 \
Stdev: 71.0996 \
5th percentile: 278 \
25th percentile: 387 \
50th percentile: 456 \
75th percentile: 476 \
95th percentile: 496 \
Mode: 474

### Stats for model_1:

Expected score: 497.976

Mean: 333.445 \
Stdev: 92.1199 \
5th percentile: 176 \
25th percentile: 267 \
50th percentile: 334 \
75th percentile: 401 \
95th percentile: 478 \
Mode: 380

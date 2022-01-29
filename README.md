# Despirala-Engine

Finds the optimal strategy for the dice game Despirala.
Can also evaluate luck and mistakes made during play.

## Rules of the game

1. At the start of each turn you roll all six dice for free and then gain five goods.
2. After that you choose to either reroll the dice (for the cost of one good) or choose some combination that you have not attempted before.
3. After picking a combination you have to complete it.
4. To that end, you repeatedly take the dice not found in the combination and reroll only them (for the cost of one good).
5. If you successfully complete the combination, you win the amount of points it is worth and your turn ends.
6. If you run out of goods before completing the chosen combination, your turn ends and you do not get any points for that combination.
7. The six single number combinations (ones, twos, ...) are an exception: you choose when you are done rolling and get points proportional to the number of correct dice.
8. At the end of the game you get one bonus point per good you have left over.

Under normal play, the goal is to maximize your score. Under misère play, the goal is to minimize your score, but rerolling is not allowed.

## List of combinations

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

## User interface

The UI is entirely in the console. It should be fairly intuitive. Moves are written in the format: name of the move, followed by its arguments given as numbers. E.g. "Collect 3", "Two triples 3 5", "Despirala", "Reroll". While collecting the two possible moves are "Continue collecting" and "Stop collecting" (or just "Continue" and "Stop"). Move names are case insensitive.

## Performance stats

Stats for normal_exact.model and misere_exact.model are obtained with 100 million tests.
All othere stats are obtained with a million tests.
Note that for an exact model, "expected score" will be equal to the empirical mean.

### Stats for normal_exact.model

Expected score: 443.616

Mean: 443.622 \
Stdev: 61.454 \
5th percentile: 310 \
25th percentile: 420 \
50th percentile: 468 \
75th percentile: 484 \
95th percentile: 501 \
Mode: 478

### Stats for normal_1000.model

Expected score: 443.694

Mean: 443.505 \
Stdev: 61.519 \
5th percentile: 310 \
25th percentile: 420 \
50th percentile: 468 \
75th percentile: 483 \
95th percentile: 501 \
Mode: 478

### Stats for normal_50.model

Expected score: 444.524

Mean: 442.461 \
Stdev: 62.156 \
5th percentile: 308 \
25th percentile: 418 \
50th percentile: 467 \
75th percentile: 483 \
95th percentile: 500 \
Mode: 477

### Stats for normal_5.model

Expected score: 457.534

Mean: 425.650 \
Stdev: 71.143 \
5th percentile: 278 \
25th percentile: 386 \
50th percentile: 456 \
75th percentile: 476 \
95th percentile: 496 \
Mode: 473

### Stats for normal_1.model

Expected score: 497.976

Mean: 334.253 \
Stdev: 91.886 \
5th percentile: 177 \
25th percentile: 268 \
50th percentile: 336 \
75th percentile: 402 \
95th percentile: 478 \
Mode: 378

### Stats for misere_exact.model

Expected score: 105.975

Mean: 105.973 \
Stdev: 52.456 \
5th percentile: 28 \
25th percentile: 71 \
50th percentile: 99 \
75th percentile: 137 \
95th percentile: 201 \
Mode: 80

### Stats for misere_50.model

Expected score: 105.857

Mean: 106.105 \
Stdev: 52.483 \
5th percentile: 28 \
25th percentile: 71 \
50th percentile: 99 \
75th percentile: 137 \
95th percentile: 201 \
Mode: 78

### Stats for misere_5.model

Expected score: 103.404

Mean: 107.960 \
Stdev: 52.653 \
5th percentile: 29 \
25th percentile: 72 \
50th percentile: 102 \
75th percentile: 139 \
95th percentile: 203 \
Mode: 78

### Stats for misere_1.model

Expected score: 99.6856

Mean: 114.095 \
Stdev: 53.207 \
5th percentile: 34 \
25th percentile: 78 \
50th percentile: 108 \
75th percentile: 146 \
95th percentile: 210 \
Mode: 85

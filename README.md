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

The UI is entirely in the console. It should be fairly intuitive. Moves are written in the format: name of the move, followed by its arguments given as numbers. E.g. "Collect 3", "Two triples 3 5", "Despirala", "Reroll". While collecting the two possible moves are "Continue collecting" and "Stop collecting" (or just "Continue" and "Stop"). Move names are case insensitive. To list all currently possible moves use "List options" (or just "List" or "Options").

# Engine

The engine precomputes the expected scores of reachable positions. If wanted, it can do this with Monte Carlo sampling in each state.
However, exactly computing the value of each state is tractable. Thus, the engine can optimally solve the game. Using these precomputed
values it can evaluate the expected score from any states and the delta in any transition between states (i.e. luck when rolling and
mistakes when choosing a move). Note that the engine solves for maximizing the expected score, without caring about other players
or winning chances.

# Competitive play

There is a WIP competitive strategy. It tries to optimize the winning chances/expected rank. It explores the possible moves of the
current player and assume the game is played according to the EV maximizing strategy from there on. It then uses the distribution
of the scores of the other players and the distributions of all the moves to choose the move which results in the best expected rank.
This is still WIP and works through Monte Carlo simulations of the games (instead of exact computations). To beat the EV maximizing
strategy it needs to simulate about 10k games per distribution, which is quite slow. Its mean rank then is around 1.443 (it wins 55.7%
of the games). This is from 4635 tests and the (one-sided) p-value is < 3e-15 (null hypothesis being that it wins 50% of the time).

## Performance stats

All stats are obtained with a million tests.
Note that for an exact model, "expected score" will actually be equal to the expected value of the score,
i.e. empirical mean will approach the calculated "expected score".

### Stats for normal_exact.model

Expected score: 443.616

Mean: 443.754 \
Stdev: 61.373 \
5th percentile: 310 \
25th percentile: 420 \
50th percentile: 468 \
75th percentile: 484 \
95th percentile: 501 \
Mode: 478

### Stats for normal_1000.model

Expected score: 443.665

Mean: 443.667 \
Stdev: 61.388 \
5th percentile: 311 \
25th percentile: 421 \
50th percentile: 468 \
75th percentile: 484 \
95th percentile: 501 \
Mode: 479

### Stats for normal_50.model

Expected score: 445.344

Mean: 442.396 \  
Stdev: 62.388 \
5th percentile: 307 \
25th percentile: 417 \
50th percentile: 467 \
75th percentile: 483 \
95th percentile: 501 \
Mode: 476

### Stats for normal_5.model

Expected score: 454.471

Mean: 429.211 \
Stdev: 70.139 \
5th percentile: 282 \
25th percentile: 392 \
50th percentile: 459 \
75th percentile: 479 \
95th percentile: 498 \
Mode: 473

### Stats for normal_1.model

Expected score: 497.152

Mean: 329.579 \
Stdev: 92.806 \
5th percentile: 179 \
25th percentile: 261 \
50th percentile: 326 \
75th percentile: 398 \
95th percentile: 477 \
Mode: 468

### Stats for misere_exact.model

Expected score: 105.973

Mean: 105.917 \
Stdev: 52.408 \
5th percentile: 28 \
25th percentile: 71 \
50th percentile: 99 \
75th percentile: 137 \
95th percentile: 201 \
Mode: 80

### Stats for misere_50.model

Expected score: 105.956

Mean: 105.889 \
Stdev: 52.443 \
5th percentile: 28 \
25th percentile: 71 \
50th percentile: 99 \
75th percentile: 137 \
95th percentile: 201 \
Mode: 80

### Stats for misere_5.model

Expected score: 103.943

Mean: 107.710 \
Stdev: 52.700 \
5th percentile: 29 \
25th percentile: 72 \
50th percentile: 101 \
75th percentile: 138 \
95th percentile: 204 \
Mode: 80

### Stats for misere_1.model

Expected score: 100.375

Mean: 113.820 \
Stdev: 53.628 \
5th percentile: 33 \
25th percentile: 77 \
50th percentile: 109 \
75th percentile: 145 \
95th percentile: 211 \
Mode: 83

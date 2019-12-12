# Despirala-Solver

Finds the optimal strategy for the dice game Despirala.

Rules of the game:

On each turn, you roll 6 dice. After seeing the result, you chose to either reroll the dice or choose some combination. If you choose a combination, it should be one you have not done before. After choosing a combination you have to complete it. To that end, you are alowed to take some of the dice (of your choosing) and reroll only them. If you successfully complete the combination you win the amount of points it is worth and your turn ends. However, you have a limitted number of rolls, called goods, and each roll consumes one good. At the start of each turn you gain 6 goods (though you consume 1 of them immediatelly for the first throw). If you run out of goods before completing the combination, your turn ends and you do not get any points for that combination. A slight exception is the six single number combinations (ones, twos, ...). With them you can choose when you have completed the combination and you get points proportional to the number of correct dice. At the end of the game you get one bouns point per good you have left over.

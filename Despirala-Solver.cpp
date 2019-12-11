#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdlib.h>
#include <time.h>

const int NUM_DICE = 6;
const int NUM_SIDES = 6;
const int MAX_ROLL = 46656; // NUM_DICE ^ NUM_SIDES
const int MAX_CODE = 1 << (NUM_SIDES + NUM_DICE);
const int DIFF_CODES = 30;
const int NUM_COMBOS = 14;
const int GOODS_PER_TURN = 5;
const double POINTS_PER_GOOD = 1;
const int MAX_GOODS = NUM_COMBOS * GOODS_PER_TURN;
const int NUM_TRIALS_1 = 100000;
const int NUM_TRIALS_2 = 10;
const int NUM_TEST_TRIALS = 100;

typedef std::array<int, NUM_SIDES> Occurs;

struct Target
{
    const std::string name;
    const std::vector<int> args;
    const int points;
    const Occurs occurs;

    Target(const std::string& name, const std::vector<int>& args, int points, const Occurs& occurs):
        name(name),
        args(args),
        points(points),
        occurs(occurs) {}
};

struct Combo
{
    Combo(const std::string& name, int points):
        name(name),
        points(points) {}
    std::string getName() const
    {
        return name;
    }
    virtual std::vector<Target> getTargets(const Occurs& diceOccurs)
    {
        return {};
    }
    virtual int getNumber() const
    {
        return 0;
    }

protected:
    const std::string name;
    const int points;
};

struct SingleCombo : Combo
{
    SingleCombo(const std::string& name, int number):
        Combo(name, 0),
        number(number) {}
    int getNumber() const
    {
        return number;
    }

protected:
    int number;
};

struct FixedCombo : Combo
{
public:
    FixedCombo(const std::string& name, int points, const Occurs& occurs):
        Combo(name, points),
        occurs(occurs) {}
    std::vector<Target> getTargets(const Occurs& diceOccurs)
    {
        return {Target(name, {}, points, occurs)};
    }

protected:
    const Occurs occurs;
};

struct PermCombo : Combo
{
public:
    PermCombo(const std::string& name, int points, const Occurs& occurs):
        Combo(name, points),
        templateOccurs(occurs)
    {
        std::sort(templateOccurs.begin(), templateOccurs.end());
    }
    std::vector<Target> getTargets(const Occurs& diceOccurs)
    {
        int code = 0;
        for (int i = 0; i < NUM_SIDES; ++i)
        {
            code = code * (NUM_DICE + 1) + diceOccurs[i];
        }
        if (cache.find(code) == cache.end())
        {
            cache[code] = rawGetTargets(diceOccurs);
        }
        return cache[code];
    }

protected:
    Occurs templateOccurs;
    std::unordered_map<int, std::vector<Target>> cache;

    virtual std::vector<Target> rawGetTargets(const Occurs& diceOccurs) const
    {
        std::array<std::pair<int, int>, NUM_SIDES> diceOccursPairs;
        makePairs(diceOccurs, diceOccursPairs);
        std::sort(diceOccursPairs.begin(), diceOccursPairs.end());

        Occurs occurs;
        std::vector<int> args = makeOccurs(templateOccurs, diceOccursPairs, occurs);

        return {Target(name, args, points, occurs)};
    }
    static void makePairs(const Occurs& diceOccurs, std::array<std::pair<int, int>, NUM_SIDES>& diceOccursPairs)
    {
        for (int i = 0; i < NUM_SIDES; ++i)
        {
            diceOccursPairs[i] = {diceOccurs[i], i};
        }
    }
    static std::vector<int> makeOccurs(const Occurs& templateOccurs, const std::array<std::pair<int, int>, NUM_SIDES>& diceOccursPairs, Occurs& currOccurs)
    {
        std::vector<int> args;
        for (int i = NUM_SIDES - 1; i >= 0; --i)
        {
            currOccurs[diceOccursPairs[i].second] = templateOccurs[i];
            if (templateOccurs[i] > 0) args.push_back(diceOccursPairs[i].second + 1);
        }
        return args;
    }
};

struct SPermCombo : PermCombo
{
    SPermCombo(const std::string& name, const Occurs& occurs):
        PermCombo(name, 0, occurs) {};

protected:
    std::vector<Target> rawGetTargets(const Occurs& diceOccurs) const
    {
        std::array<std::pair<int, int>, NUM_SIDES> diceOccursPairs;
        makePairs(diceOccurs, diceOccursPairs);
        std::sort(diceOccursPairs.begin(), diceOccursPairs.end());

        Occurs currTemplate(templateOccurs);

        int maxPoints = -1;
        std::vector<Target> ts;
        Occurs occurs;
        do
        {
            std::vector<int> args = makeOccurs(currTemplate, diceOccursPairs, occurs);
            int currPoints = evalOccurs(occurs);
            if (currPoints > maxPoints)
            {
                maxPoints = currPoints;
                ts.push_back(Target(name, args, currPoints, occurs));
            }
        }
        while (std::next_permutation(currTemplate.begin(), currTemplate.end()));

        return ts;
    }
    static int evalOccurs(const Occurs& occurs)
    {
        int points = 0;
        for (int i = 0;i < NUM_SIDES; ++i)
        {
            points += (i + 1) * occurs[i];
        }
        return points;
    }
};

Combo* combos[]={new SingleCombo("Ones",   1),
                 new SingleCombo("Twos",   2),
                 new SingleCombo("Threes", 3),
                 new SingleCombo("Fours",  4),
                 new SingleCombo("Fives",  5),
                 new SingleCombo("Sixes",  6),
                 new  SPermCombo("Three Pairs",        {2, 2, 2, 0, 0, 0}),
                 new  SPermCombo("Two Triples",        {3, 3, 0, 0, 0, 0}),
                 new   PermCombo("Four of a kind", 40, {4, 0, 0, 0, 0, 0}),
                 new  FixedCombo("Kamerun",        45, {0, 0, 0, 1, 2, 3}),
                 new  FixedCombo("Straight",       50, {1, 1, 1, 1, 1, 1}),
                 new   PermCombo("Six of a kind",  60, {6, 0, 0, 0, 0, 0}),
                 new  FixedCombo("General",        70, {0, 0, 0, 0, 0, 6}),
                 new  FixedCombo("Despirala",      80, {5, 0, 0, 0, 0, 1})};

int occursToCode(Occurs occurs)
{
    std::sort(occurs.begin(), occurs.end(), std::greater<int>());
    int code = 1;
    int curr = 0;
    while (curr < NUM_SIDES)
    {
        while (occurs[curr])
        {
            --occurs[curr];
            code *= 2;
        }
        code = code * 2 + 1;
        ++curr;
    }
    code /= 2;
    return code;
}

bool codeToOccurs(int code, Occurs& occurs)
{
    occurs.fill(0);
    int n = NUM_SIDES + NUM_DICE;
    int instrNum = 0;
    int instrs[n + 1];
    while (code && instrNum <= n)
    {
        instrs[instrNum] = code % 2;
        code /= 2;
        ++instrNum;
    }
    if (instrNum == 0 || instrNum > n) return false;
    --instrNum;
    int curr = 0;
    while (instrNum && curr < NUM_SIDES)
    {
        --instrNum;
        if (instrs[instrNum]) ++curr;
        else ++occurs[curr];
    }
    if (curr != NUM_SIDES - 1) return false;
    for (int i = 1; i < NUM_SIDES; ++i)
    {
        if (occurs[i - 1] < occurs[i]) return false;
    }
    return true;
}

void rollToOccurs(int roll, Occurs& occurs)
{
    occurs.fill(0);
    while (roll)
    {
        ++occurs[roll % NUM_SIDES];
        roll /= NUM_SIDES;
    }
}

std::unordered_map<int, int> codeToIndex;
int indexToCode[DIFF_CODES];

void generateLefts()
{
    Occurs occurs;
    int lastIdx = 0;
    for (int i = 0; i <= MAX_CODE; ++i)
    {
        if (codeToOccurs(i, occurs))
        {
            indexToCode[lastIdx] = i;
            codeToIndex[i] = lastIdx;
            ++lastIdx;
        }
    }
}

double rollsDistr[DIFF_CODES][MAX_GOODS + 1];

void findRollsDistrSingle(int idx)
{
    Occurs occurs;
    Occurs occurs2;
    int code = indexToCode[idx];
    codeToOccurs(code, occurs);
    for (int i = 0; i < NUM_TRIALS_1; ++i)
    {
        occurs2 = occurs;
        bool changed = false;
        int rolls = 0;
        while (!changed)
        {
            for (int i = 0; i < NUM_DICE; ++i)
            {
                int side = rand() % NUM_SIDES;
                if (occurs2[side])
                {
                    --occurs2[side];
                    changed = true;
                }
            }
            ++rolls;
        }
        int newCode = occursToCode(occurs2);
        int newIdx = codeToIndex[newCode];
        for (int i = 0; i <= MAX_GOODS - rolls; ++i)
        {
            rollsDistr[idx][i + rolls] += rollsDistr[newIdx][i] / NUM_TRIALS_1;
        }
    }
}

void findRollsDistr()
{
    rollsDistr[0][0] = 1;
    for (int i = 1; i < DIFF_CODES; ++i)
    {
        findRollsDistrSingle(i);
    }
}

double leftDistr[NUM_DICE + 1][NUM_DICE + 1];

int fact(int n)
{
    return n ? n * fact(n - 1) : 1;
}

void findLeftDistr()
{
    for (int i = 0; i <= NUM_DICE; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            leftDistr[i][j] = pow(5, j) / pow(6, i) * fact(i) / (fact(j) * fact(i - j));
        }
    }
}

bool isFree(int free, int idx)
{
    return (free >> idx) & 1;
}

int setUsed(int free, int idx)
{
    return free - (1 << idx);
}

bool isFound[1 << NUM_COMBOS][MAX_GOODS + 1];
double score[1 << NUM_COMBOS][MAX_GOODS + 1];
bool isFoundColl[1 << NUM_COMBOS][MAX_GOODS + 1][NUM_SIDES][NUM_DICE];
double collScore[1 << NUM_COMBOS][MAX_GOODS + 1][NUM_SIDES][NUM_DICE];

struct Move
{
    double score;
    std::string name;

    Move(double score, const std::string& name):
        score(score),
        name(name) {}
};

bool operator<(const Move& a, const Move& b)
{
    return a.score < b.score;
}

double getScore(int free, int goods);
double getScoreCont(int free, int goods);
double getCollScore(int free, int goods, int num, int left);

Move getCollMove(int free, int goods, int num, int left)
{
    Move best = Move(getScore(free, goods), "Stop");
    if (goods && left)
    {
        double sc = 0;
        for (int i = 0 ; i <= left; ++i)
        {
            sc += leftDistr[left][i] * ((left - i) * num + getCollScore(free, goods - 1, num, i));
        }
        Move option = Move(sc, "Continue");
        best = std::max(best, option);
    }
    return best;
}

double getCollScore(int free, int goods, int num, int left)
{
    if (!isFoundColl[free][goods][num][left])
    {
        isFoundColl[free][goods][num][left] = true;
        collScore[free][goods][num][left] = getCollMove(free, goods, num, left).score;
    }
    return collScore[free][goods][num][left];
}

double simulateMove(int free, int goods, const Occurs& occurs, double reward)
{
    double score = 0;
    int code = occursToCode(occurs);
    int idx = codeToIndex[code];

    double succP = 0;
    for (int i = 0; i <= goods; ++i)
    {
        double currP = rollsDistr[idx][i];
        if (currP > 0)
        {
            succP += currP;
            score += currP * getScore(free, goods - i);
        }
    }
    score += succP * reward;
    score += (1 - succP) * getScore(free, 0);
    return score;
}

Move getMove(int free, int goods, const Occurs& occurs)
{
    Move best = Move(-1, "Error: no move found");
    for (int i = 0; i < NUM_COMBOS; ++i)
    {
        if (!isFree(free, i)) continue;
        int newFree = setUsed(free, i);
        int num = combos[i]->getNumber();
        if (num)
        {
            int curr = occurs[num - 1];
            Move option = Move(num * curr + getCollScore(newFree, goods, num, curr), combos[i]->getName());
            best = std::max(best, option);
        }
        else
        {
            std::vector<Target> ts = combos[i]->getTargets(occurs);
            for (int j = 0; j < ts.size(); ++j)
            {
                const Target& t = ts[j];
                Occurs newOccurs = t.occurs;
                for (int i = 0; i < NUM_SIDES; ++i)
                {
                    newOccurs[i] = std::max(newOccurs[i] - occurs[i], 0);
                }
                Move option = Move(simulateMove(newFree, goods, newOccurs, t.points), t.name);
                best = std::max(best, option);
            }
        }
    }
    if (goods > 0)
    {
        Move option = Move(getScoreCont(free, goods - 1), "Reroll");
        best = std::max(best, option);
    }
    return best;
}

int cnt = 0;

double getScoreCont(int free, int goods)
{
    if (isFound[free][goods]) return score[free][goods];
    isFound[free][goods] = true;
    if (free == 0) score[free][goods] = goods * POINTS_PER_GOOD;
    else
    {
        Occurs occurs;
        for (int i = 0; i < NUM_TRIALS_2; ++i)
        {
            occurs.fill(0);
            for (int j = 0; j < NUM_DICE; ++j)
            {
                ++occurs[rand() % NUM_SIDES];
            }
            //rollToOccurs(rand() % MAX_ROLL, occurs);
            Move mov = getMove(free, goods, occurs);
            score[free][goods] += mov.score / NUM_TRIALS_2;
        }
    }
    ++cnt;
    //std::cerr << "score[" << free << "][" << goods << "] = " << score[free][goods] << "; // " << cnt << "\n";
    if (cnt % 10000 == 0) std::cerr << cnt << "\n";
    return score[free][goods];
}

double getScore(int free, int goods)
{
    return free ? getScoreCont(free, goods + GOODS_PER_TURN) : getScoreCont(free, goods);
}

void findExpectedScores()
{
    getScore((1 << NUM_COMBOS) - 1, 0);
}

int main()
{
    srand(time(0));
    std::cerr << "Starting." << std::endl;
    generateLefts();
    std::cerr << "Generated possible things left to roll." << std::endl;
    findRollsDistr();
    std::cerr << "Found the distribution of dice rolls needed for them." << std::endl;
    findLeftDistr();
    std::cerr << "Found the distribution of number of dice left after roll." << std::endl;
    findExpectedScores();
    std::cerr << "Found the expected scores." << std::endl;
    std::cerr << "Expected score for the game: " << getScore((1 << NUM_COMBOS) - 1, 0) << "." << std::endl;
    std::cerr << "Considered " << cnt << " cases." << std::endl;

    return 0;
}

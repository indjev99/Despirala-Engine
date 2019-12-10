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
const int MAX_GOODS = NUM_COMBOS * GOODS_PER_TURN;
const int NUM_TRIALS_1 = 100000;
const int NUM_TRIALS_2 = 5;

struct Target
{
    std::string name;
    int points;
    int occurs[NUM_SIDES];

    Target(const std::string& name, int points, const int occurs[NUM_SIDES]):
        name(name),
        points(points)
    {
        for (int i = 0; i < NUM_SIDES; ++i)
        {
            this->occurs[i] = occurs[i];
        }
    }
};

struct Combo
{
    std::string name;
    int points;

    Combo(const std::string& name, int points):
        name(name),
        points(points) {}
    virtual std::vector<Target> getTargets(const int diceOccurs[NUM_SIDES]) const = 0;
};

struct EmptyCombo : Combo
{
    EmptyCombo(const std::string& name):
        Combo(name, 0) {}
    std::vector<Target> getTargets(const int diceOccurs[NUM_SIDES]) const
    {
        return {};
    }
};

struct FixedCombo : Combo
{
    int occurs[NUM_SIDES];
    FixedCombo(const std::string& name, int points, const int occurs[NUM_SIDES]):
        Combo(name, points)
    {
        for (int i = 0; i < NUM_SIDES; ++i)
        {
            this->occurs[i] = occurs[i];
        }
    }
    std::vector<Target> getTargets(const int diceOccurs[NUM_SIDES]) const
    {
        return {Target(name, points, occurs)};
    }
};

struct PermCombo : Combo
{
    bool fixedPoints;
    int templateOccurs[NUM_SIDES];

    PermCombo(const std::string& name, int points, const int occurs[NUM_SIDES]):
        Combo(name, points)
    {
        fixedPoints = points > 0;
        for (int i = 0; i < NUM_SIDES; ++i)
        {
            templateOccurs[i] = occurs[i];
        }
        std::sort(templateOccurs, templateOccurs + NUM_SIDES);
    }
    std::vector<Target> getTargets(const int diceOccurs[NUM_SIDES])const
    {
        std::string currName = name;
        int currPoints = points;
        int occurs[NUM_SIDES];
        std::pair<int, int> diceOccurs2[NUM_SIDES];
        for (int i = 0; i < NUM_SIDES; ++i)
        {
            diceOccurs2[i] = {diceOccurs[i], i};
        }
        std::sort(diceOccurs2, diceOccurs2 + NUM_SIDES);
        for (int i = NUM_SIDES - 1; i >= 0; --i)
        {
            occurs[diceOccurs2[i].second] = templateOccurs[i];
            if (templateOccurs[i] > 0) currName += " " + std::to_string(diceOccurs2[i].second + 1);
        }
        if (!fixedPoints)
        {
            currPoints = 0;
            for (int i = 0;i < NUM_SIDES; ++i)
            {
                currPoints += (i + 1) * occurs[i];
            }
        }
        /// TO-DO: Change this so it returns all potentially better possibilities
        return {Target(currName, currPoints, occurs)};
    }
};

Combo* combos[]={new EmptyCombo("Ones"),
                 new EmptyCombo("Twos"),
                 new EmptyCombo("Threes"),
                 new EmptyCombo("Fours"),
                 new EmptyCombo("Fives"),
                 new EmptyCombo("Sixes"),
                 new  PermCombo("Three Pairs",    0,  (int[]){2, 2, 2, 0, 0, 0}),
                 new  PermCombo("Two Triples",    0,  (int[]){3, 3, 0, 0, 0, 0}),
                 new  PermCombo("Four of a kind", 40, (int[]){4, 0, 0, 0, 0, 0}),
                 new FixedCombo("Kamerun",        45, (int[]){0, 0, 0, 1, 2, 3}),
                 new FixedCombo("Straight",       50, (int[]){1, 1, 1, 1, 1, 1}),
                 new  PermCombo("Six of a kind",  60, (int[]){6, 0, 0, 0, 0, 0}),
                 new FixedCombo("General",        70, (int[]){0, 0, 0, 0, 0, 6}),
                 new FixedCombo("Despirala",      80, (int[]){5, 0, 0, 0, 0, 1})};


int occursToCode(const int occurs[NUM_SIDES])
{
    int occurs2[NUM_SIDES];
    for (int i = 0; i < NUM_SIDES; ++i)
    {
        occurs2[i] = occurs[i];
    }
    std::sort(occurs2, occurs2 + NUM_SIDES, std::greater<int>());
    int code = 1;
    int curr = 0;
    while (curr < NUM_SIDES)
    {
        while (occurs2[curr])
        {
            --occurs2[curr];
            code *= 2;
        }
        code = code * 2 + 1;
        ++curr;
    }
    code /= 2;
    return code;
}

bool codeToOccurs(int code, int occurs[NUM_SIDES])
{
    for (int i = 0; i < NUM_SIDES; ++i)
    {
        occurs[i] = 0;
    }
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

void rollToOccurs(int roll, int occurs[NUM_SIDES])
{
    for (int i = 0; i < NUM_SIDES; ++i)
    {
        occurs[i] = 0;
    }
    while (roll)
    {
        --occurs[roll % NUM_SIDES];
        roll /= NUM_SIDES;
    }
}

std::unordered_map<int, int> codeToIndex;
int indexToCode[DIFF_CODES];

void generateLefts()
{
    int occurs[NUM_SIDES];
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
    int occurs[NUM_SIDES];
    int occurs2[NUM_SIDES];
    int code = indexToCode[idx];
    codeToOccurs(code, occurs);
    for (int i = 0; i < NUM_TRIALS_1; ++i)
    {
        for (int i = 0; i < NUM_SIDES; ++i)
        {
            occurs2[i] = occurs[i];
        }
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

int newGoods(int free)
{
    return free ? GOODS_PER_TURN : 0;
}

double simulateMove(int free, int goods, const int occurs[NUM_SIDES], double reward)
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
            score += currP * getScore(free, goods - i + newGoods(free));
        }
    }
    score += succP * reward;
    score += (1 - succP) * getScore(free, newGoods(free));
    return score;
}

Move getMove(int free, int goods, const int occurs[NUM_SIDES])
{
    Move best = Move(-1, "Error: no move found");
    for (int i = 0; i < NUM_COMBOS; ++i)
    {
        if (!isFree(free, i)) continue;
        int newFree = setUsed(free, i);
        if (i < NUM_SIDES)
        {
            Move option = Move((i + 1) * occurs[i] + getScore(newFree, goods + newGoods(newFree)), combos[i]->name);
            best = std::max(best, option);
        }
        else
        {
            std::vector<Target> ts = combos[i]->getTargets(occurs);
            for (int j = 0; j < ts.size(); ++j)
            {
                Target& t = ts[j];
                for (int i = 0; i < NUM_SIDES; ++i)
                {
                    t.occurs[i] = std::max(t.occurs[i] - occurs[i], 0);
                }
                Move option = Move(simulateMove(newFree, goods, t.occurs, t.points), t.name);
                best = std::max(best, option);
            }
        }
    }
    if (goods > 0)
    {
        Move option = Move(getScore(free, goods - 1), "Reroll");
        best = std::max(best, option);
    }
    return best;
}

long long cnt = 0;

double getScore(int free, int goods)
{
    if (isFound[free][goods]) return score[free][goods];
    isFound[free][goods] = true;
    if (free == 0) score[free][goods] = goods;
    else
    {
        int occurs[NUM_SIDES];
        for (int i = 0; i < NUM_TRIALS_2; ++i)
        {
            for (int j = 0; j < NUM_SIDES; ++j)
            {
                occurs[j] = 0;
            }
            for (int j = 0; j < NUM_DICE; ++j)
            {
                ++occurs[rand() % NUM_SIDES];
            }
            Move mov = getMove(free, goods, occurs);
            score[free][goods] += mov.score / NUM_TRIALS_2;
        }
    }
    ++cnt;
    //std::cerr << "score[" << free << "][" << goods << "] = " << score[free][goods] << "; // " << cnt << "\n";
    //if (cnt % 10000 == 0) std::cerr << cnt << "\n";
    return score[free][goods];
}

void findExpectedScores()
{
    getScore((1 << NUM_COMBOS) - 1, GOODS_PER_TURN);
}

void test()
{
    int occurs[6] = {0, 2, 0, 2, 2, 0};
    const Combo* currCombo = combos[7];
    std::vector<Target> ts = currCombo->getTargets(occurs);
    for (int i = 0; i < ts.size(); i++)
    {
        std::cerr << ts[i].name << " for " << ts[i].points << " :";
        for (int j =0; j < NUM_SIDES; ++j)
        {
            std::cerr << " " << ts[i].occurs[j];
        }
        std::cerr << std::endl;
    }
}

int main()
{
    //test();

    srand(time(0));
    std::cerr << "Starting." << std::endl;
    generateLefts();
    std::cerr << "Generated possible things left to roll." << std::endl;
    findRollsDistr();
    std::cerr << "Found the distribution of dice rolls needed for them." << std::endl;
    findExpectedScores();
    std::cerr << "Found the expected scores." << std::endl;
    std::cerr << "Expected score for the game: " << score[(1 << NUM_COMBOS) - 1][GOODS_PER_TURN] << "." << std::endl;
    std::cerr << "Considered " << cnt << " cases." << std::endl;

    return 0;
}

#include <iostream>
#include <algorithm>
#include <string>
#include <stdlib.h>

const int NUM_DICE = 6;
const int NUM_SIDES = 6;
const int MAX_ROLL = 46656; // NUM_DICE ^ NUM_SIDES
const int MAX_CODE = 1 << (NUM_SIDES + NUM_DICE);
const int DIFF_CODES = 30;
const int MAX_GOODS = 60;
const int NUM_TRIALS = 100000;

struct Combo
{
    std::string name;
    int points;
    int occurs[NUM_SIDES];

    Combo(const std::string& name, int points, const int occurs[NUM_SIDES]):
        name(name),
        points(points)
    {
        for (int i = 0; i < NUM_SIDES; ++i)
        {
            this->occurs[i] = occurs[i];
        }
    }
    virtual void set_dice(const int diceOccurs[NUM_SIDES], int goods) {}
};

struct PermCombo : Combo
{
    bool fixedPoints;
    int templateOccurs[NUM_SIDES];
    PermCombo(const std::string& name, int points, const int occurs[NUM_SIDES]):
        Combo(name, points, occurs)
    {
        fixedPoints = points > 0;
        for (int i = 0; i < NUM_SIDES; ++i)
        {
            templateOccurs[i] = occurs[i];
        }
        std::sort(templateOccurs, templateOccurs + NUM_SIDES);
    }
    void set_dice(const int diceOccurs[NUM_SIDES], int goods)
    {
        std::pair<int, int> diceOccurs2[NUM_SIDES];
        for (int i = 0; i < NUM_SIDES; ++i)
        {
            diceOccurs2[i] = {diceOccurs[i], i};
        }
        std::sort(diceOccurs2, diceOccurs2 + NUM_SIDES);
        for (int i = 0; i < NUM_SIDES; ++i)
        {
            occurs[i] = templateOccurs[diceOccurs2[i].second];
        }
        if (!fixedPoints)
        {
            points = 0;
            for (int i = 0;i < NUM_SIDES; ++i)
            {
                points += (i + 1) * occurs[i];
            }
        }
    }
};

Combo* combos[]={new PermCombo("Three Pairs",    0,  (int[]){2, 2, 2, 0, 0, 0}),
                 new PermCombo("Two Triples",    0,  (int[]){3, 3, 0, 0, 0, 0}),
                 new PermCombo("Four of a kind", 40, (int[]){4, 0, 0, 0, 0, 0}),
                 new     Combo("Kamerun",        45, (int[]){0, 0, 0, 1, 2, 3}),
                 new     Combo("Straight",       50, (int[]){1, 1, 1, 1, 1, 1}),
                 new PermCombo("Six of a kind",  60, (int[]){6, 0, 0, 0, 0, 0}),
                 new     Combo("General",        70, (int[]){0, 0, 0, 0, 0, 6}),
                 new     Combo("Despirala",      80, (int[]){5, 0, 0, 0, 0, 1})};


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

int codeToIndex[MAX_CODE + 1];
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
    for (int i = 0; i < NUM_TRIALS; ++i)
    {
        for (int i = 0; i < NUM_SIDES; ++i)
        {
            occurs2[i] = occurs[i];
        }
        int newCode = code;
        int rolls = 0;
        while (newCode == code)
        {
            for (int i = 0; i < NUM_DICE; ++i)
            {
                int side = rand() % 6;
                if (occurs2[side]) -- occurs2[side];
            }
            ++rolls;
            newCode = occursToCode(occurs2);
        }
        int newIdx = codeToIndex[newCode];
        for (int i = 0; i <= MAX_GOODS - rolls; ++i)
        {
            rollsDistr[idx][i + rolls] += rollsDistr[newIdx][i] / NUM_TRIALS;
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

int main()
{
    srand(0);
    generateLefts();
    std::cerr << "Generated possible things left to roll." << std::endl;
    findRollsDistr();
    std::cerr << "Found the distribution of dice rolls needed for them." << std::endl;


    return 0;
}

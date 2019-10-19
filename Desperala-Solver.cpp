#include <iostream>
#include <algorithm>
#include <string>

const int NUM_DICE=6;
const int NUM_SIDES=6;
const int MAX_ROLL=46656; //6^6
const int DIFF_ROLLS=462; //11! / (6! * 5!)

struct Roll
{
    int occurs[NUM_SIDES]; //occur
    int code;
    int index;

    Roll() {}

    //code - between 0 and 6^6
    Roll(int code)
    {
        for (int i=0;i<NUM_SIDES;++i)
        {
            occurs[i]=0;
        }
        int dice[NUM_DICE];
        for (int i=0;i<NUM_DICE;++i)
        {
            dice[i]=code%NUM_SIDES;
            code/=NUM_SIDES;
            ++occurs[dice[i]];
        }
        std::sort(dice, dice+NUM_DICE);
        code=0;
        for (int i=0;i<NUM_DICE;++i)
        {
            code*=NUM_SIDES;
            code+=dice[i];
        }
        this->code=code;
    }
};

Roll rolls[MAX_ROLL];
int indexToCode[DIFF_ROLLS];
int codeToIndex[MAX_ROLL];

void generate_rolls()
{
    int sortedCodes[MAX_ROLL];
    for (int i=0;i<MAX_ROLL;++i)
    {
        rolls[i]=Roll(i);
        sortedCodes[i]=rolls[i].code;
    }
    std::sort(sortedCodes, sortedCodes+MAX_ROLL);
    int lastIdx=-1;
    int lastCode=-1;
    for (int i=0;i<MAX_ROLL;++i)
    {
        if (sortedCodes[i]>lastCode)
        {
            lastCode=sortedCodes[i];
            ++lastIdx;
            indexToCode[lastIdx]=lastCode;
            codeToIndex[lastCode]=lastIdx;
        }
    }
}

struct Combo
{
    std::string name;
    int points;

    Combo(const std::string& name, int points, const int occurs[NUM_SIDES]):
        name(name),
        points(points)
    {
        for (int i=0;i<NUM_SIDES;++i)
        {
            this->occurs[i]=occurs[i];
        }
    }
    virtual void set_dice(const int diceOccurs[NUM_SIDES]) {}
};

struct FixedCombo : Combo
{
    int occurs[NUM_SIDES];
    FixedCombo(const std::string& name, int points, const int occurs[NUM_SIDES]):
        name(name),
        points(points)
    {
        for (int i=0;i<NUM_SIDES;++i)
        {
            this->occurs[i]=occurs[i];
        }
    }
    virtual void set_dice(const int diceOccurs[NUM_SIDES]) {}
};

struct PermCombo : Combo
{
    bool fixedPoints;
    int templateOccurs[NUM_SIDES];
    PermCombo(const std::string& name, int points, const int occurs[NUM_SIDES]):
        Combo(name, points, occurs)
    {
        fixedPoints=points>0;
        for (int i=0;i<NUM_SIDES;++i)
        {
            templateOccurs[i]=occurs[i];
        }
        std::sort(templateOccurs, templateOccurs+NUM_SIDES);
    }
    void set_dice(const int diceOccurs[NUM_SIDES])
    {
        std::pair<int, int> diceOccurs2[NUM_SIDES];
        for (int i=0;i<NUM_SIDES;++i)
        {
            diceOccurs2[i]={diceOccurs[i], i};
        }
        std::sort(diceOccurs2, diceOccurs2+NUM_SIDES);
        for (int i=0;i<NUM_SIDES;++i)
        {
            occurs[i]=templateOccurs[diceOccurs2[i].second];
        }
        if (!fixedPoints)
        {
            points=0;
            for (int i=0;i<NUM_SIDES;++i)
            {
                points+=(i+1)*occurs[i];
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
                 new     Combo("Desperala",      80, (int[]){5, 0, 0, 0, 0, 1})};

int main()
{
    generate_rolls();

    return 0;
}

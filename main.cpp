#include <iostream>
#include "programs.h"
#include <time.h>

using namespace std;


int main()
{
    srand(time(NULL));
    int number;
    cout<<"1 -\ttime Phrase plus Laypunov\n2 -\tbiffur and Laypunov\n3-\ttime Phrase plus Laypunov sd count\n\n";
    cin>>number;
    switch (number)
    {
        case 1:
        {
            timePhraseLaypunov();
        }
        break;
        case 2:
        {
            biffurAndLambda();
        }
        break;
        case 3:
        {
            ;
        }
        break;
    }

    return 0;
}

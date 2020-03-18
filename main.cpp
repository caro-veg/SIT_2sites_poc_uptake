#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <random>
#include <cmath>
#include <chrono>
#include <vector>

//Windows Issue - mingwin doesn't recognize to_string(const T& n) function
namespace patch
{
    template<typename T> std::string to_string(const T& n)
    {
        std::ostringstream stm;
        stm << n;
        return stm.str();
    }
}

using namespace std;


int main(int argc, char* argv[])
{
    cout << "Hello world!" << endl;

    //get current time to calculate run time after program ends
    chrono::time_point<chrono::system_clock> startTime, endTime, startTimeI, endTimeI;
    startTime = chrono::system_clock::now();

    //seed for MT
    unsigned seed = 12345;
    //generate a Mersenne Twister with seed
    mt19937 mtgenerator(seed);
    //generate uniform real distribution
    uniform_real_distribution<double> unifDistribution(0.0, 1.0);

    double beta, f, gamma, g1, g2, sigmaBase1, sigmaBase2, sigmaTreat1, sigmaTreat2;
    double sensitivity, specificity, u;
    int S, S_start, I00, I00_start, I10, I10_start, I01, I01_start, I11, I11_start, T00, T10, T01, T11;
    int Talt00, Talt10, Talt01, Talt11;

    double tmax; //simulation time step
    int numberOfRepeats; //number of simulation runs and number or time steps

    //read in parameter values and initial variable values from command line
    beta = strtod(argv[1], NULL);
    f = strtod(argv[2], NULL);
    gamma = strtod(argv[3], NULL);
    g1 = strtod(argv[4], NULL);
    g2 = strtod(argv[5], NULL);
    sensitivity = strtod(argv[6], NULL);
    specificity = strtod(argv[7], NULL);
    u = strtod(argv[8], NULL);
    sigmaBase1 = strtod(argv[9], NULL);
    sigmaBase2 = strtod(argv[10], NULL);
    sigmaTreat1 = strtod(argv[11], NULL);
    sigmaTreat2 = strtod(argv[12], NULL);
    S_start = strtol(argv[13], NULL, 10);
    I00_start = strtol(argv[14], NULL, 10);
    I10_start = strtol(argv[15], NULL, 10);
    I01_start = strtol(argv[16], NULL, 10);
    I11_start = strtol(argv[17], NULL, 10);
    //nobody is treated at beginning so all Tij=0

    tmax = strtod(argv[18], NULL);
    numberOfRepeats = strtol(argv[19], NULL, 10);

    cout << "beta: " << beta << "\t" << "f: " << f << "\t" << "gamma: " << gamma << "\t" << "g1: " << g1 << "\t" << "g2: " << g2 << endl;
    cout << "sensitivity: " << sensitivity << "\t" << "specificity: " << specificity << "\t" << "uptake: " << u << endl;
    cout << "sigmaBase1: " << sigmaBase1 << "\t" << "sigmaBase2: " << sigmaBase2 << "\t" << "sigmaTreat1: " << sigmaTreat1 << "\t" <<  "sigmaTreat2: " << sigmaTreat2 << endl;
    cout << "S_start: " << S_start << "\tI00_start: " << I00_start << "\tI10_start: " << I10_start << "\tI01_start: " << I01_start << "\tI11_start: " << I11_start << endl << endl;

    int N = S_start + I00_start + I10_start + I01_start + I11_start;
    //cout << "N: " << N << endl << endl;

    ofstream outfile;

    outfile.open(("SIT_2Sites-" + patch::to_string(beta) + "-" + patch::to_string(f) + "-" + patch::to_string(gamma) + "-" + patch::to_string(g1) + "-" + patch::to_string(sensitivity) + "-" + patch::to_string(specificity) + "-" + patch::to_string(u) + "-" + patch::to_string(sigmaBase1) + "-" + patch::to_string(sigmaBase2) + "-" + patch::to_string(sigmaTreat1) + "-" + patch::to_string(sigmaTreat2) + "-" + patch::to_string(S_start) + "-" + patch::to_string(I00_start) + "-" + patch::to_string(I10_start) + "-" + patch::to_string(I01_start) + "-" + patch::to_string(I11_start) + "-" + patch::to_string(tmax) + "-" + patch::to_string(numberOfRepeats) + ".csv").c_str());
    outfile << "repeat,time,S,I00,I10,I01,I11,T00,T10,T01,T11,Ta00,Ta10,Ta01,Ta11" << endl;

    /*int independentResistanceI = 0;
    int independentResistanceT = 0;
    int independentSingleI = 0;
    int independentSingleT = 0;*/

    //start simulation
    for(int i=0; i<numberOfRepeats; ++i)
    {
        startTimeI = chrono::system_clock::now();

        cout << endl << i << endl << endl;
        //set variables to initial values
        S = S_start;
        I00 = I00_start;
        I10 = I10_start;
        I01 = I01_start;
        I11 = I11_start;
        T00 = 0;
        T10 = 0;
        T01 = 0;
        T11 = 0;
        Talt00 = 0;
        Talt10 = 0;
        Talt01 = 0;
        Talt11 = 0;

        //cout << "S: " << S << "\tI00: " << I00 << "\tI10: " << I10 << "\tI01: " << I01 << "\tI11: " << I11 << endl;
        //cout << "T00: " << T00 << "\tT10: " << T10 << "\tT01: " << T01 << "\tT11: " << T11 << endl << endl;

        double t = 0;   //current time
        double tn = 0;  //next time point at which the state of the system will be recorded

        /*int independentMutantsI10 = 0;
        int independentMutantsI01 = 0;
        int independentMutantsI11 = 0;
        int independentMutantsT10 = 0;
        int independentMutantsT01 = 0;
        int independentMutantsT11 = 0;*/

        //int transgression = 0;

        do
        {
            //cout << tn << endl;
            //if it is time update the state of the system
            if(t>=tn)
            {
                outfile << i << "," << t << "," << S << "," << I00 << "," << I10 << "," << I01 << "," << I11 << "," << T00 << "," << T10 << "," << T01 << "," << T11 << "," << Talt00 << "," << Talt10 << "," << Talt01 << "," << Talt11 << endl;
                tn = tn + 0.1;
            }

            //determine the total rate of each event
            double infectionI00 = beta*S*I00;
            //cout << "infectionI00: " << infectionI00 << endl;
            double infectionI10 = beta*S*I10;
            //cout << "infectionI10: " << infectionI10 << endl;
            double infectionI01 = beta*S*I01;
            //cout << "infectionI01: " << infectionI01 << endl;
            double infectionI11 = beta*S*I11;
            //cout << "infectionI11: " << infectionI11 << endl;

            double recoveryI00 = f*I00;
            //cout << "recoveryI00: " << recoveryI00 << endl;
            double recoveryI10 = f*I10;
            //cout << "recoveryI10: " << recoveryI10 << endl;
            double recoveryI01 = f*I01;
            //cout << "recoveryI01: " << recoveryI01 << endl;
            double recoveryI11 = f*I11;
            //cout << "recoveryI11: " << recoveryI11 << endl;

            double treatmentI00 = (u*sensitivity + (1-u))*gamma*I00;
            //cout << "treatmentI00: " << treatmentI00 << endl;
            double treatmentI10 = (u*(1 - specificity) + (1-u))*gamma*I10;
            //cout << "treatmentI10: " << treatmentI10 << endl;
            double treatmentI01 = (u*(1 - specificity) + (1-u))*gamma*I01;
            //cout << "treatmentI01: " << treatmentI01 << endl;
            double treatmentI11 = (u*(1 - specificity) + (1-u))*gamma*I11;
            //cout << "treatmentI11: " << treatmentI11 << endl;

            double treatmentAltI00 = u*(1 - sensitivity)*gamma*I00;
            double treatmentAltI10 = u*specificity*gamma*I10;
            double treatmentAltI01 = u*specificity*gamma*I01;
            double treatmentAltI11 = u*specificity*gamma*I11;

            double cureT00 = g1*T00;
            //cout << "cureT00: " << cureT00 << endl;
            double cureT10 = g1*T10;
            //cout << "cureT10: " << cureT10 << endl;
            double cureT01 = g1*T01;
            //cout << "cureT01: " << cureT01 << endl;
            double failToCureT11 = g1*T11;
            //cout << "failToCureT11: " << failToCureT11 << endl;

            double cureTalt00 = g2*Talt00;
            double cureTalt10 = g2*Talt10;
            double cureTalt01 = g2*Talt01;
            double cureTalt11 = g2*Talt11;

            double mutationI00toI10 = sigmaBase1*I00;
            //cout << "mutationI00toI10: " << mutationI00toI10 << endl;
            double mutationI00toI01 = sigmaBase2*I00;
            //cout << "mutationI00toI01: " << mutationI00toI01 << endl;
            double mutationI10toI11 = sigmaBase2*I10;
            //cout << "mutationI10toI11: " << mutationI10toI11 << endl;
            double mutationI01toI11 = sigmaBase1*I01;
            //cout << "mutationI01toI11: " << mutationI01toI11 << endl;

            double mutationI10toI00 = sigmaBase1*I10;
            //cout << "mutationI10toI00: " << mutationI10toI11 << endl;
            double mutationI01toI00 = sigmaBase2*I01;
            //cout << "mutationI01toI00: " << mutationI01toI11 << endl;
            double mutationI11toI10 = sigmaBase2*I11;
            //cout << "mutationI11toI10: " << mutationI11toI10 << endl;
            double mutationI11toI01 = sigmaBase1*I11;
            //cout << "mutationI11toI01: " << mutationI11toI01 << endl;

            double mutationT00toT10 = sigmaTreat1*T00;
            //cout << "mutationT00toT10: " << mutationT00toT10 << endl;
            double mutationT00toT01 = sigmaTreat2*T00;
            //cout << "mutationT00toT01: " << mutationT00toT01 << endl;
            double mutationT10toT11 = sigmaTreat2*T10;
            //cout << "mutationT10toT11: " << mutationT10toT11 << endl;
            double mutationT01toT11 = sigmaTreat1*T01;
            //cout << "mutationT01toT11: " << mutationT01toT11 << endl << endl;


            //Determine the total event rate
            double totalRate = infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI01 + recoveryI10 + recoveryI11
                                + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI00 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                                + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10 + cureTalt01 + cureTalt11
                                + mutationI00toI10 + mutationI00toI01 + mutationI10toI11 + mutationI01toI11 + mutationI10toI00 + mutationI01toI00 + mutationI11toI10 + mutationI11toI01
                                + mutationT00toT10 + mutationT00toT01 + mutationT10toT11 + mutationT01toT11;

            //cout << "totalRate: " << totalRate << endl << endl;

            //Determine the time step at which next event occurs and update time
            //draw random double from uniform distribution [0,1]
            double rand1 = unifDistribution(mtgenerator);
            //Use random number to determine next time step
            //Exception - cannot divide by 0 or take log 0
            try
            {
                if(totalRate==0 || rand1==0)
                    throw 0;
                t += -1/totalRate * log(rand1);
            }
            catch(int e)
            {
                cout << "cannot divide by 0 or take log of 0" << endl;
                break;
            }

            //Determine the next event to occur and update the populations
            //draw random double from uniform distribution [0,1]
            double rand2 = unifDistribution(mtgenerator);
            //calculate fraction of total rate with random number
            double q = rand2 * totalRate;

            //Determine next event according to probability
            if(q <= infectionI00)
            {
                S--;
                if(S<0)
                {
                    S=0;
                    //transgression++;
                }

                I00++;
                if(I00>N)
                {
                    I00=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10))
            {
                S--;
                if(S<0)
                {
                    S=0;
                    //transgression++;
                }

                I10++;
                if(I10>N)
                {
                    I10=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01))
            {
                S--;
                if(S<0)
                {
                    S=0;
                   //transgression++;
                }
                I01++;
                if(I01>N)
                {
                    I01=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11))
            {
                S--;
                if(S<0)
                {
                    S=0;
                    //transgression++;
                }
                I11++;
                if(I11>N)
                {
                    I11=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00))
            {
                S++;
                if(S>N)
                {
                    S=N;
                    //transgression++;
                }
                I00--;
                if(I00<0)
                {
                    I00=0;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10))
            {
                S++;
                if(S>N)
                {
                    S=N;
                    //transgression++;
                }
                I10--;
                if(I10<0)
                {
                    I10=0;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01))
            {
                S++;
                if(S>N)
                {
                    S=N;
                    //transgression++;
                }
                I01--;
                if(I01<0)
                {
                    I01=0;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11))
            {
                S++;
                if(S>N)
                {
                    S=N;
                    //transgression++;
                }
                I11--;
                if(I11<0)
                {
                    I11=0;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00))
            {
                I00--;
                if(I00<0)
                {
                    I00=0;
                    //transgression++;
                }
                T00++;
                if(T00>N)
                {
                    T00=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10))
            {
                I10--;
                if(I10<0)
                {
                    I10=0;
                    //transgression++;
                }
                T10++;
                if(T10>N)
                {
                    T10=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01))
            {
                I01--;
                if(I01<0)
                {
                    I01=0;
                    //transgression++;
                }
                T01++;
                if(T01>N)
                {
                    T01=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11))
            {
                I11--;
                if(I11<0)
                {
                    I11=0;
                    //transgression++;
                }
                T11++;
                if(T11>N)
                {
                    T11=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00))
            {
                I00--;
                if(I00<0)
                {
                    I00=0;
                    //transgression++;
                }
                Talt00++;
                if(Talt00>N)
                {
                    Talt00=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10))
            {
                I10--;
                if(I10<0)
                {
                    I10=0;
                    //transgression++;
                }
                Talt10++;
                if(Talt10>N)
                {
                    Talt10=N;
                    //transgression++;
                }
            }
             else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01))
            {
                I01--;
                if(I01<0)
                {
                    I01=0;
                    //transgression++;
                }
                Talt01++;
                if(Talt01>N)
                {
                    Talt01=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11))
            {
                I11--;
                if(I11<0)
                {
                    I11=0;
                    //transgression++;
                }
                Talt11++;
                if(Talt11>N)
                {
                    Talt11=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00))
            {
                T00--;
                if(T00<0)
                {
                    T00=0;
                    //transgression++;
                }
                S++;
                if(S>N)
                {
                    S=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10))
            {
                T10--;
                if(T10<0)
                {
                    T10=0;
                    //transgression++;
                }
                S++;
                if(S>N)
                {
                    S=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01))
            {
                T01--;
                if(T01<0)
                {
                    T01=0;
                    //transgression++;
                }
                S++;
                if(S>N)
                {
                    S=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11))
            {
                T11--;
                if(T11<0)
                {
                    T11=0;
                    //transgression++;
                }
                I11++;
                if(I11>N)
                {
                    I11=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00))
            {
                Talt00--;
                if(Talt00<0)
                {
                    Talt00=0;
                    //transgression++;
                }
                S++;
                if(S>N)
                {
                    S=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10))
            {
                Talt10--;
                if(Talt10<0)
                {
                    Talt10=0;
                    //transgression++;
                }
                S++;
                if(S>N)
                {
                    S=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10 + cureTalt01))
            {
                Talt01--;
                if(Talt01<0)
                {
                    Talt01=0;
                    //transgression++;
                }
                S++;
                if(S>N)
                {
                    S=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10 + cureTalt01 + cureTalt11))
            {
                Talt11--;
                if(Talt11<0)
                {
                    Talt11=0;
                    //transgression++;
                }
                S++;
                if(S>N)
                {
                    S=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10 + cureTalt01 + cureTalt11
                          + mutationI00toI10))
            {
                I00--;
                if(I00<0)
                {
                    I00=0;
                    //transgression++;
                }
                I10++;
                if(I10>N)
                {
                    I10=N;
                    //transgression++;
                }

                //independentMutantsI10++;
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10 + cureTalt01 + cureTalt11
                          + mutationI00toI10 + mutationI00toI01))
            {
                I00--;
                if(I00<0)
                {
                    I00=0;
                    //transgression++;
                }
                I01++;
                if(I01>N)
                {
                    I01=N;
                    //transgression++;
                }

                //independentMutantsI01++;
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10 + cureTalt01 + cureTalt11
                          + mutationI00toI10 + mutationI00toI01 + mutationI10toI11))
            {
                I10--;
                if(I10<0)
                {
                    I10=0;
                    //transgression++;
                }
                I11++;
                if(I11>N)
                {
                    I11=N;
                    //transgression++;
                }

                //independentMutantsI11++;
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10 + cureTalt01 + cureTalt11
                          + mutationI00toI10 + mutationI00toI01 + mutationI10toI11 + mutationI01toI11))
            {
                I01--;
                if(I01<0)
                {
                    I01=0;
                    //transgression++;
                }
                I11++;
                if(I11>N)
                {
                    I11=N;
                    //transgression++;
                }

                //independentMutantsI11++;
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10 + cureTalt01 + cureTalt11
                          + mutationI00toI10 + mutationI00toI01 + mutationI10toI11 + mutationI01toI11
                          + mutationI10toI00))
            {
                I10--;
                if(I10<0)
                {
                    I10=0;
                    //transgression++;
                }
                I00++;
                if(I00>N)
                {
                    I00=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10 + cureTalt01 + cureTalt11
                          + mutationI00toI10 + mutationI00toI01 + mutationI10toI11 + mutationI01toI11
                          + mutationI10toI00 + mutationI01toI00))
            {
                I01--;
                if(I01<0)
                {
                    I01=0;
                    //transgression++;
                }
                I00++;
                if(I00>N)
                {
                    I00=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10 + cureTalt01 + cureTalt11
                          + mutationI00toI10 + mutationI00toI01 + mutationI10toI11 + mutationI01toI11
                          + mutationI10toI00 + mutationI01toI00 + mutationI11toI10))
            {
                I11--;
                if(I11<0)
                {
                    I11=0;
                    //transgression++;
                }
                I10++;
                if(I10>N)
                {
                    I10=N;
                    //transgression++;
                }
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10 + cureTalt01 + cureTalt11
                          + mutationI00toI10 + mutationI00toI01 + mutationI10toI11 + mutationI01toI11
                          + mutationI10toI00 + mutationI01toI00 + mutationI11toI10 + mutationI11toI01))
            {
                I11--;
                if(I11<0)
                {
                    I11=0;
                    //transgression++;
                }
                I01++;
                if(I01>N)
                {
                    I01=N;
                    //transgression++;
                }

            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10 + cureTalt01 + cureTalt11
                          + mutationI00toI10 + mutationI00toI01 + mutationI10toI11 + mutationI01toI11
                          + mutationI10toI00 + mutationI01toI00 + mutationI11toI10 + mutationI11toI01
                          + mutationT00toT10))
            {
                T00--;
                if(T00<0)
                {
                    T00=0;
                    //transgression++;
                }
                T10++;
                if(T10>N)
                {
                    T10=N;
                    //transgression++;
                }

                //independentMutantsT10++;
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10 + cureTalt01 + cureTalt11
                          + mutationI00toI10 + mutationI00toI01 + mutationI10toI11 + mutationI01toI11
                          + mutationI10toI00 + mutationI01toI00 + mutationI11toI10 + mutationI11toI01
                          + mutationT00toT10 + mutationT00toT01))
            {
                T00--;
                if(T00<0)
                {
                    T00=0;
                    //transgression++;
                }
                T01++;
                if(T01>N)
                {
                    T01=N;
                    //transgression++;
                }

                //independentMutantsT01++;
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10 + cureTalt01 + cureTalt11
                          + mutationI00toI10 + mutationI00toI01 + mutationI10toI11 + mutationI01toI11
                          + mutationI10toI00 + mutationI01toI00 + mutationI11toI10 + mutationI11toI01
                          + mutationT00toT10 + mutationT00toT01 + mutationT10toT11))
            {
                T10--;
                if(T10<0)
                {
                    T10=0;
                    //transgression++;
                }
                T11++;
                if(T11>N)
                {
                    T11=N;
                    //transgression++;
                }

                //independentMutantsT11++;
            }
            else if(q <= (infectionI00 + infectionI10 + infectionI01 + infectionI11 + recoveryI00 + recoveryI10 + recoveryI01 + recoveryI11
                          + treatmentI00 + treatmentI10 + treatmentI01 + treatmentI11 + treatmentAltI00 + treatmentAltI10 + treatmentAltI01 + treatmentAltI11
                          + cureT00 + cureT10 + cureT01 + failToCureT11 + cureTalt00 + cureTalt10 + cureTalt01 + cureTalt11
                          + mutationI00toI10 + mutationI00toI01 + mutationI10toI11 + mutationI01toI11
                          + mutationI10toI00 + mutationI01toI00 + mutationI11toI10 + mutationI11toI01
                          + mutationT00toT10 + mutationT00toT01 + mutationT10toT11 + mutationT01toT11))
            {
                T01--;
                if(T01<0)
                {
                    T01=0;
                    //transgression++;
                }
                T11++;
                if(T11>N)
                {
                    T11=N;
                    //transgression++;
                }

                //independentMutantsT11++;
            }


            //terminate simulation if it has reached tmax or GC has died out
        }while(t<=tmax && ((I11 + T11 + Talt11) <= 0.1*(I00 + I10 + I01 + I11 + T00 + T10 + T01 + T11 + Talt00 + Talt10 + Talt01 + Talt11)) && ((I00 + I10 + I01 + I11 + T00 + T10 + T01 + T11 + Talt00 + Talt10 + Talt01 + Talt11)>0));

        //get current time and calculate runtime of iteration i
        endTimeI = chrono::system_clock::now();
        chrono::duration<double> elapsed_secondsI = endTimeI - startTimeI;
        //output runtime to console
        cout << "iteration " << i << ":\t" << elapsed_secondsI.count() << endl << endl;
        //cout << "transgressions:\t" << transgression << endl;
        /*cout << "I00 --> I10:\t" << independentMutantsI10 << endl;
        cout << "I00 --> I01:\t" << independentMutantsI01 << endl;
        cout << "I10/I01 --> I11:\t" << independentMutantsI11 << endl;
        cout << "T00 --> T10:\t" << independentMutantsT10 << endl;
        cout << "T00 --> T01:\t" << independentMutantsT01 << endl;
        cout << "T10/T01 --> T11:\t" << independentMutantsT11 << endl;
        cout << endl;*/

        /*independentResistanceI += independentMutantsI11;
        independentResistanceT += independentMutantsT11;
        independentSingleI = independentSingleI + independentMutantsI01 + independentMutantsI10;
        independentSingleT = independentSingleT + independentMutantsT01 + independentMutantsT10;*/
    }

    //close output file
    outfile.close();

    //get current time and calculate total runtime of programme
    endTime = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = endTime - startTime;
    //output runtime to console
    cout << elapsed_seconds.count() << endl << endl;

    /*cout << "mean independent single mutants untreated: " << static_cast<double>(independentSingleI)/numberOfRepeats << endl;
    cout << "mean independent single mutants treated: " << static_cast<double>(independentSingleT)/numberOfRepeats << endl;
    cout << "mean independent resistance mutants untreated: " << static_cast<double>(independentResistanceI)/numberOfRepeats << endl;
    cout << "mean independent resistance mutants treated: " << static_cast<double>(independentResistanceT)/numberOfRepeats << endl << endl;*/

    return 0;
}

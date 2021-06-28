#include <iostream>
#include <vector>
#include <eigen3/Eigen/Eigen>
#include <fstream>

using namespace std;
using namespace Eigen;

//printing progress during calculation
void print(auto a)
{
    std::cout << a << std::endl;
}

//display progress during simulation
void display_progress(double progress)
{
    int width = 70;
    cout << "[";
    int position = width * progress;
    for (int i = 0; i < width; ++i)
    {
        if (i < position)
            cout << "=";
        else if (i == position)
            cout << ">";
        else
            cout << " ";
    }
    cout << "] " << int(progress * 100.0) << " %\r";
    cout.flush();
}

// =================================================================================================
//                      PROGRAMMSTRUKTUR
//
//                          ===========
//                          | Dataset |
//                          ===========
//                              |
//                           ========
//                           | Data |
//                           ========
//                              |
// ==============             ======                =============
// | Thermostat | ----------- | MD | -------------- | Potential |
// ==============             ======                =============
//                              |
//                           ========
//                           | main |
//                           ========
//
// - Die Klasse MD beinhaltet die primäre Logik des Algorithmus
// - Die Klasse Thermostat und Potential sind davon getrennt, damit unterschiedliche
//   Thermostate und Potentiale durch Vererbung implementiert und flexibel verwendet
//   werden können
// - Die Klasse Data speichert die in der MD-Simulation gespeicherten Daten und
//   kümmert sich um das Abspeichern
//      - Dataset ist eine Datensatz aus Zeit, Temperatur, ...
//      - Data hält mehrere Datensätze und noch einige Daten, die nicht zeitaufgelöst
//        gespeichert werden
//      - Statt umständlich getter und setter zu verwenden, sind die Member von Data
//        und Dataset public, da es sich um simple Datencontainer handelt
// - main() ruft MD mit den für die Aufgabenteile notwendigen Parametern auf
//
// Hinweise zu den verwendeten Vectoren:
//      - Wegen der Performance verwenden wir Vector2d statt VectorXd
//        Das macht jedoch den möglichen Übergang zu 3d-Zuständen umständlicher als mit VectorXd
//      - Für Listen von Daten wird std::vector werden
// =================================================================================================

// ================================ Potential-Klasse ===============================================

// Virtuelle Klasse, aus der konkrete Potentiale vererbt werden können
// (Hier nur Lennard-Jones nötig, aber so kann man schnell andere Potentiale implementieren)
class Potential
{
public:
    virtual double V(double r2) const = 0;    // Virtuelle Funktion
    virtual Vector2d F(Vector2d r) const = 0; // Virtuelle Funktion
};

class PotentialLJ : public Potential
{
public:
    double V(double r2) const;    // Überschreibt virtuelle Funktion
    Vector2d F(Vector2d r) const; // Überschreibt virtuelle Funktion
};

// Für Potential reicht das Quadrat der Vektorlänge, womit eine Wurzelberechnung gespart wird
double PotentialLJ::V(double r2) const
{
    double r6_inv = pow(r2, -3);

    return 4 * r6_inv * (r6_inv - 1.);
}

Vector2d PotentialLJ::F(Vector2d r) const
{
    double r2_inv = 1. / r.squaredNorm();
    return 48 * r * pow(r2_inv, 4) * (pow(r2_inv, 3) - 0.5);
}

// ------------------------------ Ende Potential-Klasse ------------------------------------------
// ================================ Thermostat-Klasse ===============================================

// Virtuelle Klasse, aus der konkrete Thermostate vererbt werden können
class Thermostat
{
public:
    virtual void rescale(vector<Vector2d> &v, double T) const = 0;
};

// Kein Thermostat
class NoThermostat : public Thermostat
{
public:
    void rescale(vector<Vector2d> &v, double T) const {} // Macht nichts
};

// Isokinetisches Thermostat für Aufgabe d)
class IsokinThermostat : public Thermostat
{
public:
    void rescale(vector<Vector2d> &v, double T) const;
};

void IsokinThermostat::rescale(vector<Vector2d> &v, double T) const
{
    double tmp_sum = 0;
    for (Vector2d i : v)
    {
        tmp_sum += i.dot(i);
    }
    double alpha = sqrt(30 * T / tmp_sum);

    for (Vector2d &i : v)
    {
        i *= alpha;
    }
}

// ------------------------------ Ende Thermostat-Klasse ------------------------------------------
// ================================ Data-Structs ===============================================

// Datensatz für zeitaufgelöste Daten
// (structs im Grunde gleich zu class, aber alle Member sind standardmäßig public)
struct Dataset
{
    double t, T, Ekin, Epot;
    Vector2d vS;
    MatrixXd R; //used to create the animation
};

// Rückgabedaten der MD-Simulation
// Konstruktor Data data(n); reserviert Speicher und füllt Paarkorrelationsfunktion mit 0en
struct Data
{
    vector<Dataset> datasets; // Zeitaufgelöste Datensätze
    vector<double> rBin, g;   // Gemittelte Paarkorrelationsfunktion
    // vector<Vector2d> r;       // Momentaufnahme der finalen Position
    //    Für Aufgabe e) kann es sinnvoll sein, r stattdessen
    //    in den zeitaufgelösten Datasets abzuspeichern

    vector<double> equi_tmp; // saves the cumulative values of ekin, epot, T during the equilibration
                             // is needed, to calculate the zeitgemittelte observables afterwards

    Data(uint n, uint numBins, double binSize);
    void save(const string &filenameSets,
              const string &filenameG,
              const string &filenameR) const;
};

Data::Data(uint n, uint numBins, double binSize) : datasets(n), // Initializer list, weil sie Konstruktoren der Member aufruft
                                                   rBin(numBins),
                                                   g(numBins, 0.),
                                                   //    r(0),
                                                   equi_tmp(4)
{
}

void Data::save(const string &filenameSets, const string &filenameG, const string &filenameR) const
{
    ofstream outfile;

    //used to create the animation
    outfile.open("build/anim." + filenameSets, ios::out);
    for (Dataset i : datasets)
    {
        outfile << i.R << "\n";
    }
    outfile.close();

    outfile.open("build/" + filenameSets, ios::out);

    outfile << "#" << "t" << "\t" << "vS" << "\t" << "Ekin" << "\t" << "Epot" << "\t" << "T" << "\n";

    for (Dataset i : datasets)
    {
        outfile << i.t << "\t" << i.vS.norm() << "\t" << i.Ekin << "\t" << i.Epot << "\t" << i.T << "\n";
    }
    outfile.close();

    outfile.open("build/" + filenameG, ios::out);
    outfile << "#g"
            << "\n";
    for (double i : g)
    {
        outfile << i << "\n";
    }
    outfile.close();

    outfile.open("build/" + filenameR, ios::out);
    outfile << "#r"
            << "\n";
    for (double i : rBin)
    {
        outfile << i << "\n";
    }
    outfile.close();
}

// ------------------------------ Ende Data-Structs ------------------------------------------
// ================================ MD-Klasse ===============================================

class MD
{
public:
    MD(double L, uint N, uint particlesPerRow, double T,
       Potential &potential, Thermostat &thermostat,
       uint numBins = 1000);

    void equilibrate(const double dt, const unsigned int n);
    Data measure(const double dt, const unsigned int n);

private:
    vector<Vector2d> r, v;
    double L;
    uint N;
    Potential &potential;
    Thermostat &thermostat;
    double T;
    double t = 0.;

    //e_pot is calculated in calcAcc
    double e_pot_acc = 0;

    uint numBins;
    double binSize;

    // Teilchen werden in Box [0,L]x[0,L] verschoben
    void centerParticles();

    // Berechnungen wichtiger Messgrößen
    double calcT() const;
    double calcEkin() const;
    double calcEpot() const;

    Vector2d calcvS() const;

    Dataset calcDataset() const;

    // Berechnung der Beschleunigung
    //  Um redundante Rechnungen zu vermeiden, kann es sinnvoll sein, das Histogram
    //  bei der Berechnung der Beschleunigungen zu aktualisieren, daher wird es
    //  hier als Referenz übergeben
    vector<Vector2d> calcAcc(vector<double> &hist); // const;

    // Berechnung des Abstandsvektors zwischen Teilchen r[i] und nähesten Spiegelteilchen von r[j]
    Vector2d calcDistanceVec(uint i, uint j) const;
};

// Initialisierung des Systems per Konstruktor
MD::MD(double L, uint N, uint particlesPerRow, double T,
       Potential &potential, Thermostat &thermostat,
       uint numBins) : L(L),
                       N(N),
                       T(T),
                       potential(potential),
                       thermostat(thermostat),
                       numBins(numBins),
                       binSize(L / (2 * numBins))
{
    //initialize grid and random velocities
    for (int i = 0; i < particlesPerRow; i++)
    {
        for (int j = 0; j < particlesPerRow; j++)
        {
            Vector2d r_tmp;
            r_tmp(0) = 1. / 8. * (1. + 2 * i) * L;
            r_tmp(1) = 1. / 8. * (1. + 2 * j) * L;
            r.push_back(r_tmp);

            Vector2d v_tmp = Vector2d::Random();
            v.push_back(v_tmp);
        }
    }

    //cms velocity
    Vector2d vS = calcvS();
    for (Vector2d &i : v)
    {
        i -= vS;
    }

    //scale with given temperature, starting temperature
    double sum_tmp = 0;
    for (Vector2d i : v)
    {
        sum_tmp += i.dot(i);
    }

    double alpha = sqrt((2 * 16 - 2) * T / sum_tmp);
    for (Vector2d &i : v)
    {
        i *= alpha;
    }
}

// Integration ohne Datenaufnahme zur reinen Äquilibrierung
void MD::equilibrate(const double dt, const unsigned int n)
{
    //during eq. no need to calc the histogram. a vector of size one yields one iteration step for l in calcAcc
    vector<double> dummy(1);

    //displays progress
    double progress = 0.0;
    double progress_increment = 1 / double(n);

    vector<Vector2d> a_tmp_1 = calcAcc(dummy);

    for (int i = 0; i < n; i++)
    {
        //displays progress
        display_progress(progress);
        progress += progress_increment;

        //tmp acceleration for verlet algorithm
        vector<Vector2d> a_tmp_0 = a_tmp_1;

        for (int j = 0; j < N; j++)
        {
            r[j] += dt * v[j] + 0.5 * dt * dt * a_tmp_0[j];
        }

        centerParticles();

        a_tmp_1 = calcAcc(dummy);

        for (int j = 0; j < N; j++)
        {
            v[j] += 0.5 * (a_tmp_1[j] + a_tmp_0[j]) * dt;
        }
    }
    //display last 100% and next line
    display_progress(progress_increment * n);
    cout << "\n";
}

Data MD::measure(const double dt, const unsigned int n)
{
    Data dataframe(n, numBins, binSize);

    vector<double> histogram(numBins);

    //displays progress
    double progress = 0.0;
    double progress_increment = 1 / double(n);

    vector<Vector2d> a_tmp_1 = calcAcc(histogram);

    for (int i = 0; i < n; i++)
    {
        t++;

        //displays progress
        display_progress(progress);
        progress += progress_increment;

        vector<Vector2d> a_tmp_0 = a_tmp_1;

        for (int j = 0; j < N; j++)
        {
            r[j] += dt * v[j] + 0.5 * dt * dt * a_tmp_0[j];
        }

        centerParticles();

        //hist is calculated in accAcc
        a_tmp_1 = calcAcc(histogram);

        for (int j = 0; j < N; j++)
        {
            v[j] += 0.5 * (a_tmp_1[j] + a_tmp_0[j]) * dt;
        }
        thermostat.rescale(v, T);

        dataframe.datasets[i] = calcDataset();
    }

    //display last 100% and next line
    display_progress(progress_increment * n);
    cout << "\n";

    //calc of pair correlation function
    for (int l = 0; l < numBins; l++)
    {
        double del_V = M_PI * (pow(binSize * (l + 1), 2) - pow(binSize * l, 2));
        double rho = N / (L * L);
        dataframe.g[l] += histogram[l] / (N * rho * del_V * n);
        dataframe.rBin[l] = l * binSize;
    }

    return dataframe;
}

void MD::centerParticles()
{
    for (Vector2d &i : r)
    {
        i(0) -= L * floor(i(0) / L);
        i(1) -= L * floor(i(1) / L);
    }
}

double MD::calcT() const
{
    return calcEkin() / 15; // T = 2 Ekin / k_B N_f,  N_f = 2*N-2 = 30
}

double MD::calcEkin() const
{
    double Ekin = 0;

    for (Vector2d i : v)
    {
        Ekin += i.dot(i);
    }
    return 0.5 * Ekin;
}

double MD::calcEpot() const
{
    return e_pot_acc; //normally 1/2 Epot, but because of Newtons third law, two energies can be calculated at once each loop
}

Vector2d MD::calcvS() const
{
    Vector2d vS(0, 0);
    for (Vector2d i : v)
    {
        vS += i;
    }
    return vS / N;
}

Dataset MD::calcDataset() const
{
    Dataset tmp_data;

    tmp_data.Ekin = calcEkin();
    tmp_data.Epot = e_pot_acc;
    tmp_data.T = calcT();
    tmp_data.vS = calcvS();
    tmp_data.t = t;

    //creating animations
    MatrixXd tmp_mat(2, N);
    for (int j = 0; j < N; j++)
    {
        tmp_mat.col(j) = r[j];
    }
    tmp_data.R = tmp_mat;

    return tmp_data;
}

Vector2d MD::calcDistanceVec(uint i, uint j) const
{
    double L_half_square = L * L / 4;
    Vector2d r_ij(0, 0);

    for (int m = -1; m <= 1; m++)
    {
        for (int n = -1; n <= 1; n++)
        {
            r_ij = r[i] - (r[j] + Vector2d(m * L, n * L));

            if (r_ij.squaredNorm() < L_half_square) //use squaredNorm to save one sqrt-calc
            {
                return r_ij;
            }
        }
    }
    return r_ij; //if no particle is in Radius L/2 -> return (0, 0)
}

vector<Vector2d> MD::calcAcc(vector<double> &hist)
{
    e_pot_acc = 0;
    vector<Vector2d> a(N, Vector2d::Zero());
    double bin_size = hist.size();
    bool loop_switch = true; //because acc has to be calculated only once

    for (int l = 0; l < bin_size; l++) //for pair correlation calculation (one iteration when equilibrating)
    {
        double r_min = l * binSize;
        double r_max = (l + 1) * binSize;

        for (int i = 0; i < N - 1; i++)
        {
            for (int j = i + 1; j < N; j++)
            {
                Vector2d r_ij = calcDistanceVec(i, j);
                if (r_min <= r_ij.norm() && r_ij.norm() < r_max)
                {
                    hist[l] += 2;
                }

                if (r_ij.squaredNorm() > 1e-12 && loop_switch == true) //avoid numerical instabilities , i.e. r^2 = 1e-15 -> should be zero
                {
                    Vector2d F_tmp = potential.F(r_ij);
                    a[i] += F_tmp;
                    a[j] -= F_tmp;

                    e_pot_acc += potential.V(r_ij.squaredNorm());
                }
            }
        }
        loop_switch = false;
    }
    return a;
}

// ------------------------------ Ende MD-Klasse ------------------------------------------

int main(void)
{
    PotentialLJ LJ;
    NoThermostat noThermo;
    IsokinThermostat isoThermo;

    const uint partPerRow = 4;
    const uint N = 16;
    const double L = 8;
    const int numBins = 100;
    // b) Äquilibrierungstest
    {
        const double T = 1.; // T = 2*E_kin / k_B*N_f = 1/15
        const double dt = 0.01;
        const uint steps = 1e5;

        MD md(L, N, partPerRow, T, LJ, noThermo, numBins);

        print("b): T=1, equilibration");
        md.measure(dt, steps).save("b)set.dat", "b)g.dat", "b)r.dat");
    }

    // c) Paarkorrelationsfunktion
    {
        string TstringVec[3] = {"0.01", "1", "100"};
        vector<double> dt{0.01, 0.01, 0.001};

        for (int i = 0; i < 3; i++)
        {
            const double T = stod(TstringVec[i]);
            const uint equiSteps = 1e5;
            const uint steps = 1e5;

            MD md(L, N, partPerRow, T, LJ, noThermo, numBins);

            print("c): T=" + TstringVec[i] + ", equilibration");
            md.equilibrate(dt[i], equiSteps);

            print("c): T=" + TstringVec[i] + ", measurement");
            md.measure(dt[i], steps).save("c)set" + TstringVec[i] + ".dat", "c)g" + TstringVec[i] + ".dat", "c)r" + TstringVec[i] + ".dat");
        }
    }
    // d) Thermostat
    {
        const double T = 0.01;
        const double steps = 1e5;
        const double dt = 0.01;
        MD md(L, N, partPerRow, T, LJ, isoThermo, numBins);

        print("d): T=0.01, equilibration");
        md.measure(dt, steps).save("d)set.dat", "d)g.dat", "d)r.dat");
    }

    return 0;
}

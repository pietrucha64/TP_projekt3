#define _USE_MATH_DEFINES
#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <cmath>
#include <vector>
#include <set>
#include <AudioFile.h>


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)


namespace py = pybind11;


int sinus(double frequency)
{
    using namespace matplot;
    
    std::vector<double> x = linspace(0, 4 * M_PI);
    std::vector<double> y = transform(x, [frequency](auto x) { return sin(frequency*x); });
    plot(x, y);
    title("Sygnal Sinusoidalny");
    xlabel("x");
    ylabel("sin(frequency*x)");
    show();
    return 0;
}

int cosinus(double frequency)
{
    using namespace matplot;
    
    std::vector<double> x = linspace(0, 4 * M_PI);
    std::vector<double> y = transform(x, [frequency](auto x) { return cos(frequency*x); });
    plot(x, y);
    title("Sygnal Cosinusoidalny");
    xlabel("x");
    ylabel("cos(frequency*x)");
    show();
    return 0;
}

int prostokatny(double frequency)
{
    using namespace matplot;
    double f = 1/frequency;
    std::vector<double> x;
    std::vector<double> y;
    for(int k = 0; k < 20;k+=2)
    {
        x.push_back(f*k);
        y.push_back(-1);
        x.push_back(f*k);
        y.push_back(1);
        x.push_back(f*k+f);
        y.push_back(1);
        x.push_back(f*k+f);
        y.push_back(-1);
        x.push_back(f*k + 2*f);
        y.push_back(-1);
    }
    plot(x, y);
    title("Sygnal Prostokatny");
    xlabel("X");
    ylabel("Y");
    show();
    return 0;
}

int piloksztaltny(double frequency)
{
    using namespace matplot;
    double f = 1/frequency;
    std::vector<double> x;
    std::vector<double> y;
    for(int k = 0; k < 10;k++)
    {
        x.push_back(f*k);
        y.push_back(-1);
        x.push_back(f*k +f);
        y.push_back(1);
    }
    plot(x, y);
    title("Sygnal Piloksztaltny");
    xlabel("X");
    ylabel("Y");
    show();
    return 0;
}
//Wizualizacja sygnalu audio
int Wizualizacja()
{
    using namespace matplot;
    std::vector<double> x;
    std::vector<double> y;

    const std::string inputFilePath = std::string (PROJECT_BINARY_DIR) + "/test-audio.wav";
    AudioFile<double> audioFile;
    audioFile.load (inputFilePath);
    int channel = 0;
    int numSamples = audioFile.getNumSamplesPerChannel();
    if(numSamples>400)
    {
        numSamples = 400;
    }
        
    for (int i = 0; i < numSamples; i++)
    {
        double currentSample = audioFile.samples[channel][i];
        x.push_back(i);
        y.push_back(currentSample);
    }
    plot(x, y);
    title("Wizualizacja Sygnalu");
    xlabel("X");
    ylabel("Y");

    return 0;
}
//Progowanie sygnalu audio
int progowanie(double prog){

    using namespace matplot;
    std::vector<double> x;
    std::vector<double> y;

    const std::string inputFilePath = std::string (PROJECT_BINARY_DIR) + "/test-audio.wav";
    AudioFile<double> audioFile;
    audioFile.load (inputFilePath);
    int channel = 0;
    int numSamples = audioFile.getNumSamplesPerChannel();
    if(numSamples>400)
    {
        numSamples = 400;
    }
        
    for (int i = 0; i < numSamples; i++)
    {
        double currentSample = audioFile.samples[channel][i];
        x.push_back(i);
        if(currentSample > prog)
            y.push_back(1);
        else
            y.push_back(0);

    }
    plot(x, y);
    title("progowanie");
    xlabel("X");
    ylabel("Y");
  
    return 0;
}
//DFT dla sinus(frequency)
int DFT(double frequency)
{
    using namespace matplot;
    
    std::vector<double> x = linspace(0, 4*frequency * M_PI);
    std::vector<double> a = transform(x, [frequency](auto x) { return sin(frequency*x); });
    std::vector<double> REX;
    std::vector<double> IMX;

    for (int k = 0; k < a.size(); k++)
    {
        double A = 0;
        double B = 0;
        for(int n = 0; n < a.size(); n++)
        {
            A+= cos((2 * M_PI * k * n) / a.size()) * a[n];
            B+= -sin((2 * M_PI * k * n) / a.size()) * a[n];
        }
        REX.push_back(A);
        IMX.push_back(B);
    }
    plot(x,IMX);
    title("DFT(sin)");
    xlabel("X");
    ylabel("IMX");

    return 0;
}
//IDFT dla sin
int IDFT(double frequency)
{
    using namespace matplot;
    
    std::vector<double> x = linspace(0, 4*frequency * M_PI);
    std::vector<double> a = transform(x, [frequency](auto x) { return sin(frequency*x); });
    std::vector<double> REX;
    std::vector<double> IMX;

    for (int k = 0; k < a.size(); k++)
    {
        double A = 0;
        double B = 0;
        for(int n = 0; n < a.size(); n++)
        {
            A+= cos((2 * M_PI * k * n) / a.size()) * a[n];
            B+= -sin((2 * M_PI * k * n) / a.size()) * a[n];
        }
        REX.push_back(A);
        IMX.push_back(B);
    }
    
    std::vector<double> Y;
    for (int i = 0; i < a.size(); i++)
    {
        double c = 0;
        for (int k = 0; k < a.size(); k++)
        {
            c += REX[k] * cos(2*M_PI*k*i/a.size());
            c += -IMX[k] * sin(2*M_PI*k*i/a.size());
        }
        Y.push_back(c);
        
    }

    plot(x,Y);
    title("IDFT");
    xlabel("X");
    ylabel("Y");
    return 0;
}

PYBIND11_MODULE(TP_projekt3, m) {

    m.def("sinus",&sinus);
    m.def("cosinus",&cosinus);
    m.def("prostokatny",&prostokatny);
    m.def("piloksztaltny",&piloksztaltny);
    m.def("Wizualizacja",&Wizualizacja);
    m.def("progowanie",&progowanie);
    m.def("DFT",&DFT);
    m.def("IDFT",&IDFT);

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
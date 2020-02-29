

#include <ql/qldefines.hpp>
#ifdef BOOST_MSVC
#include <ql/auto_link.hpp>
#endif
#include <ql/instruments/vanillaoption.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/pricingengines/vanilla/analytichestonengine.hpp>
#include <ql/pricingengines/vanilla/baroneadesiwhaleyengine.hpp>
#include <ql/pricingengines/vanilla/batesengine.hpp>
#include <ql/pricingengines/vanilla/binomialengine.hpp>
#include <ql/pricingengines/vanilla/bjerksundstenslandengine.hpp>
#include <ql/pricingengines/vanilla/fdamericanengine.hpp>
#include <ql/pricingengines/vanilla/fdbermudanengine.hpp>
#include <ql/pricingengines/vanilla/fdeuropeanengine.hpp>
#include <ql/pricingengines/vanilla/integralengine.hpp>
#include <ql/pricingengines/vanilla/mcamericanengine.hpp>
#include <ql/pricingengines/vanilla/mceuropeanengine.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/utilities/dataformatters.hpp>

#include <boost/timer.hpp>
#include <iomanip>
#include <iostream>

#include <ql/experimental/lattices/extendedbinomialtree.hpp>
#include <ql/termstructures/volatility/equityfx/blackvariancecurve.hpp>
//#include <ql/quantlib.hpp>
using namespace std;
using namespace QuantLib;

int main() {

    try {
boost::timer timer;
std::cout << std::endl;
Calendar calendar = TARGET();
Date todaysDate(27, Jan, 2020);
Date settlementDate(27, Jan, 2020);
Settings::instance().evaluationDate() = todaysDate;

// option parameters
Option::Type type(Option::Call);
Real stock = 35;
Real strike = 40;
Spread dividendYield = 0.00;
Rate riskFreeRate = 0.05;
Volatility volatility = 0.20;
Date maturity(27, May, 2019);
DayCounter dayCounter = Actual365Fixed();

boost::shared_ptr<Exercise> europeanExercise(new EuropeanExercise(maturity));

Handle<Quote> underlyingH(boost::shared_ptr<Quote>(new SimpleQuote(stock)));

// bootstrap the yield/dividend/vol curves
Handle<YieldTermStructure> flatTermStructure(boost::shared_ptr<YieldTermStructure>(
	new FlatForward(
	settlementDate,
	riskFreeRate,
	dayCounter)));

Handle<YieldTermStructure> flatDividendTS(boost::shared_ptr<YieldTermStructure>(
	new FlatForward(settlementDate,
	dividendYield,
	dayCounter)));

Handle<BlackVolTermStructure> flatVolTS(boost::shared_ptr<BlackVolTermStructure>(
	new BlackConstantVol(
	settlementDate,
	calendar,
	volatility,
	dayCounter)));

Handle<ImpliedVolTermStructure> impVolTS(boosImpliedVolTermStructure(constRelinkableHandle<BlackVolTermStructure>&originalTS, constDate&newReferenceDate)

boost::shared_ptr<StrikedTypePayoff> payoff(
	new PlainVanillaPayoff(
	type,
	strike));

boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(
	new BlackScholesMertonProcess(
	underlyingH,
	flatDividendTS,
	flatTermStructure,
	flatVolTS));


// our option is European-style
VanillaOption europeanOption(
	payoff,
	europeanExercise);

// computing the option price with the analytic Black-Scholes formulae
//europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
//	new AnalyticEuropeanEngine(
//	bsmProcess)));

// outputting
std::cout << "Option type = " << type << std::endl;
std::cout << "Maturity = " << maturity << std::endl;
std::cout << "Stock price = " << stock << std::endl;
std::cout << "Strike = " << strike << std::endl;
//std::cout << "Risk-free interest rate = " << io::rate(riskFreeRate) << std::endl;
//std::cout << "Dividend yield = " << io::rate(dividendYield) << std::endl;
//std::cout << "Volatility = " << io::volatility(volatility) << std::endl << std::endl;
//std::cout<<"European Option value = " << europeanOption.NPV() << std::endl;
Size widths[] = {45, 14, 14, 14};
    Size timeSteps = 10000;
    timer.restart();
    double seconds = timer.elapsed();
    std::cout << " \nStart in ";
    std::cout << seconds << " s\n" << std::endl;


    // Binomial method: Jarrow-Rudd
    europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
        new BinomialVanillaEngine<ExtendedJarrowRudd>(bsmProcess, timeSteps)));
    std::cout << std::setw(widths[1]) << std::left << europeanOption.NPV()
              << std::endl;

    seconds = timer.elapsed();
    std::cout << " \nRun completed in ";
    std::cout << seconds << " s\n" << std::endl;


        //std::cout << "Risk-free interest rate = " << io::rate(riskFreeRate)<< std::endl;
        //std::cout << "Dividend yield = " << io::rate(dividendYield)<< std::endl;
        //std::cout << "Volatility = " << io::volatility(volatility)<< std::endl;



//         bootstrap the yield/dividend/vol curves
//
//        Handle<YieldTermStructure> flatTermStructure(
//
//            ext::shared_ptr<YieldTermStructure>(
//
//                new FlatForward(settlementDate, riskFreeRate, dayCounter)));
//
//        Handle<YieldTermStructure> flatDividendTS(
//
//            ext::shared_ptr<YieldTermStructure>(
//
//                new FlatForward(settlementDate, dividendYield, dayCounter)));
//
//        Handle<BlackVolTermStructure> flatVolTS(
//
//            ext::shared_ptr<BlackVolTermStructure>(
//
//                new BlackConstantVol(settlementDate, calendar, volatility,
//
//                                     dayCounter)));
//
//        ext::shared_ptr<StrikedTypePayoff> payoff(
//
//                                        new PlainVanillaPayoff(type, strike));
//
//        ext::shared_ptr<BlackScholesMertonProcess> bsmProcess(
//
//                 new BlackScholesMertonProcess(underlyingH, flatDividendTS,
//
//                                               flatTermStructure, flatVolTS));



        return 0;

    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "unknown error" << std::endl;
        return 1;
    }
}


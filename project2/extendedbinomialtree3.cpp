/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2003 Ferdinando Ametrano
 Copyright (C) 2005 StatPro Italia srl
 Copyright (C) 2008 John Maiden
 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/
 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.
 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

// #include <ql/experimental/lattices/extendedbinomialtree.hpp>
#include "extendedbinomialtree3.hpp"
#include <ql/math/distributions/binomialdistribution.hpp>

namespace QuantLib {

    ExtendedJarrowRudd3::ExtendedJarrowRudd3(
                        const ext::shared_ptr<StochasticProcess1D>& process,
                        Time end, Size steps, Real)
    : ExtendedEqualProbabilitiesBinomialTree3<ExtendedJarrowRudd3>(
                                                        process, end, steps) {
        // drift removed
        for (Size i = 0; i <= steps; i ++) {
            upStepStorage.push_back(this->upStep(i));
        }
        up_ = process->stdDeviation(0.0, x0_, dt_);
    }

    Real ExtendedJarrowRudd3::upStep(Size i) const {
        Time stepTime = i*this->dt_;
        return this->treeProcess_->stdDeviation(stepTime, x0_, dt_);
    }



    ExtendedCoxRossRubinstein3::ExtendedCoxRossRubinstein3(
                        const ext::shared_ptr<StochasticProcess1D>& process,
                        Time end, Size steps, Real)
    : ExtendedEqualJumpsBinomialTree3<ExtendedCoxRossRubinstein3>(
                                                        process, end, steps) {
        for (Size i = 0; i <= steps; i ++) {
            dxStepStorage.push_back(this->dxStep(i));
            probUpStorage.push_back(this->probUp(i));
        }
        dx_ = process->stdDeviation(0.0, x0_, dt_);
        pu_ = 0.5 + 0.5*this->driftStepStorage[0] / dx_;
        pd_ = 1.0 - pu_;
        QL_REQUIRE(pu_<=1.0, "negative probability");
        QL_REQUIRE(pu_>=0.0, "negative probability");
    }

    Real ExtendedCoxRossRubinstein3::dxStep(Size i) const {
        Time stepTime = i*this->dt_;
        return this->treeProcess_->stdDeviation(stepTime, x0_, dt_);
    }

    Real ExtendedCoxRossRubinstein3::probUp(Size i) const {
        return 0.5 + 0.5*this->driftStepStorage[i]/this->dxStepStorage[i];
    }


    ExtendedAdditiveEQPBinomialTree3::ExtendedAdditiveEQPBinomialTree3(
                        const ext::shared_ptr<StochasticProcess1D>& process,
                        Time end, Size steps, Real)
    : ExtendedEqualProbabilitiesBinomialTree3<ExtendedAdditiveEQPBinomialTree3>(
                                                        process, end, steps) {
          Real driftStep_ = this->driftStepStorage[0];
          for (Size i = 0; i <= steps; i ++) {
              upStepStorage.push_back(this->upStep(i));
          }
          up_ = -0.5 * driftStep_ + 0.5 *
				      std::sqrt(4.0*process->variance(0.0, x0_, dt_) -
				          3.0*driftStep_*driftStep_);
    }

    Real ExtendedAdditiveEQPBinomialTree3::upStep(Size i) const {
      Time stepTime = i*this->dt_;
      Real driftStep_ = this->driftStepStorage[i];
      return (-0.5 * driftStep_ + 0.5 *
        std::sqrt(4.0*this->treeProcess_->variance(stepTime, x0_, dt_) -
          3.0*driftStep_*driftStep_));
    }




    ExtendedTrigeorgis3::ExtendedTrigeorgis3(
                        const ext::shared_ptr<StochasticProcess1D>& process,
                        Time end, Size steps, Real)
    : ExtendedEqualJumpsBinomialTree3<ExtendedTrigeorgis3>(process, end, steps) {
        Real driftStep_ = this->driftStepStorage[0];
        for (Size i = 0; i <= steps; i ++) {
            dxStepStorage.push_back(this->dxStep(i));
            probUpStorage.push_back(this->probUp(i));
        }
        dx_ = std::sqrt(process->variance(0.0, x0_, dt_) +
            driftStep_*driftStep_);
        pu_ = 0.5 + 0.5*driftStep_ / this->dxStepStorage[0];
        pd_ = 1.0 - pu_;

        QL_REQUIRE(pu_<=1.0, "negative probability");
        QL_REQUIRE(pu_>=0.0, "negative probability");
    }

    Real ExtendedTrigeorgis3::dxStep(Size i) const {
        Time stepTime = i*this->dt_;
        Real driftStep_ = this->driftStepStorage[i];
        return std::sqrt(this->treeProcess_->variance(stepTime, x0_, dt_) +
			     driftStep_*driftStep_);
    }

    Real ExtendedTrigeorgis3::probUp(Size i) const {
        return 0.5 + 0.5*this->driftStepStorage[i]/this->dxStepStorage[i];
    }

    ExtendedTian3::ExtendedTian3(
                        const ext::shared_ptr<StochasticProcess1D>& process,
                        Time end, Size steps, Real)
    : ExtendedBinomialTree3<ExtendedTian3>(process, end, steps) {

        Real q = std::exp(process->variance(0.0, x0_, dt_));
        Real r = std::exp(this->driftStepStorage[0])*std::sqrt(q);

        up_ = 0.5 * r * q * (q + 1 + std::sqrt(q * q + 2 * q - 3));
        down_ = 0.5 * r * q * (q + 1 - std::sqrt(q * q + 2 * q - 3));

        pu_ = (r - down_) / (up_ - down_);
        pd_ = 1.0 - pu_;


        QL_REQUIRE(pu_<=1.0, "negative probability");
        QL_REQUIRE(pu_>=0.0, "negative probability");
    }

    Real ExtendedTian3::underlying(Size i, Size index) const {
        Time stepTime = i*this->dt_;
        Real q = std::exp(this->treeProcess_->variance(stepTime, x0_, dt_));
        Real r = std::exp(this->driftStepStorage[i])*std::sqrt(q);

        Real up = 0.5 * r * q * (q + 1 + std::sqrt(q * q + 2 * q - 3));
        Real down = 0.5 * r * q * (q + 1 - std::sqrt(q * q + 2 * q - 3));

        return x0_ * std::pow(down, Real(BigInteger(i)-BigInteger(index)))
            * std::pow(up, Real(index));
    }

    Real ExtendedTian3::probability(Size i, Size, Size branch) const {
        Time stepTime = i*this->dt_;
        Real q = std::exp(this->treeProcess_->variance(stepTime, x0_, dt_));
        Real r = std::exp(this->driftStepStorage[i])*std::sqrt(q);
        Real up = 0.5 * r * q * (q + 1 + std::sqrt(q * q + 2 * q - 3));
        Real down = 0.5 * r * q * (q + 1 - std::sqrt(q * q + 2 * q - 3));

        Real pu = (r - down) / (up - down);
        Real pd = 1.0 - pu;

        return (branch == 1 ? pu : pd);
    }





    ExtendedLeisenReimer3::ExtendedLeisenReimer3(
                        const ext::shared_ptr<StochasticProcess1D>& process,
                        Time end, Size steps, Real strike)
    : ExtendedBinomialTree3<ExtendedLeisenReimer3>(process, end,
                                                 (steps%2 ? steps : steps+1)),
      end_(end), oddSteps_(steps%2 ? steps : steps+1), strike_(strike) {

        QL_REQUIRE(strike>0.0, "strike " << strike << "must be positive");
        Real driftStep_ = this->driftStepStorage[0];
        Real variance = process->variance(0.0, x0_, end);

        Real ermqdt = std::exp(driftStep_ + 0.5*variance / oddSteps_);
        Real d2 = (std::log(x0_ / strike) + driftStep_*oddSteps_) /
            std::sqrt(variance);

        pu_ = PeizerPrattMethod2Inversion(d2, oddSteps_);
        pd_ = 1.0 - pu_;
        Real pdash = PeizerPrattMethod2Inversion(d2+std::sqrt(variance),
                                                 oddSteps_);
        up_ = ermqdt * pdash / pu_;
        down_ = (ermqdt - pu_ * up_) / (1.0 - pu_);

    }

    Real ExtendedLeisenReimer3::underlying(Size i, Size index) const {
        Time stepTime = i*this->dt_;
        Real variance = this->treeProcess_->variance(stepTime, x0_, end_);

        Real driftStep_ = this->driftStepStorage[i];

        Real ermqdt = std::exp(driftStep_ + 0.5*variance / oddSteps_);
		    Real d2 = (std::log(x0_ / strike_) + driftStep_*oddSteps_) /
			     std::sqrt(variance);

        Real pu = PeizerPrattMethod2Inversion(d2, oddSteps_);
        Real pdash = PeizerPrattMethod2Inversion(d2+std::sqrt(variance),
            oddSteps_);
        Real up = ermqdt * pdash / pu;
        Real down = (ermqdt - pu * up) / (1.0 - pu);

        return x0_ * std::pow(down, Real(BigInteger(i)-BigInteger(index)))
            * std::pow(up, Real(index));
    }

    Real ExtendedLeisenReimer3::probability(Size i, Size, Size branch) const {
        Time stepTime = i*this->dt_;
        Real variance = this->treeProcess_->variance(stepTime, x0_, end_);
        Real d2 = (std::log(x0_ / strike_) + this->driftStepStorage[i]*oddSteps_) /
			std::sqrt(variance);

        Real pu = PeizerPrattMethod2Inversion(d2, oddSteps_);
        Real pd = 1.0 - pu;

        return (branch == 1 ? pu : pd);
    }



    Real ExtendedJoshi43::computeUpProb(Real k, Real dj) const {
        Real alpha = dj/(std::sqrt(8.0));
        Real alpha2 = alpha*alpha;
        Real alpha3 = alpha*alpha2;
        Real alpha5 = alpha3*alpha2;
        Real alpha7 = alpha5*alpha2;
        Real beta = -0.375*alpha-alpha3;
        Real gamma = (5.0/6.0)*alpha5 + (13.0/12.0)*alpha3
            +(25.0/128.0)*alpha;
        Real delta = -0.1025 *alpha- 0.9285 *alpha3
            -1.43 *alpha5 -0.5 *alpha7;
        Real p =0.5;
        Real rootk= std::sqrt(k);
        p+= alpha/rootk;
        p+= beta /(k*rootk);
        p+= gamma/(k*k*rootk);
        p+= delta/(k*k*k*rootk);
        return p;
    }

    ExtendedJoshi43::ExtendedJoshi43(
                        const ext::shared_ptr<StochasticProcess1D>& process,
                        Time end, Size steps, Real strike)
    : ExtendedBinomialTree3<ExtendedJoshi43>(process, end,
                                           (steps%2 ? steps : steps+1)),
      end_(end), oddSteps_(steps%2 ? steps : steps+1), strike_(strike) {

        QL_REQUIRE(strike>0.0, "strike " << strike << "must be positive");
        Real variance = process->variance(0.0, x0_, end);
        Real driftStep_ = this->driftStepStorage[0];
        Real ermqdt = std::exp(driftStep_ + 0.5*variance / oddSteps_);
		Real d2 = (std::log(x0_ / strike) + driftStep_*oddSteps_) /
			std::sqrt(variance);

        pu_ = computeUpProb((oddSteps_-1.0)/2.0,d2 );
        pd_ = 1.0 - pu_;
        Real pdash = computeUpProb((oddSteps_-1.0)/2.0,d2+std::sqrt(variance));
        up_ = ermqdt * pdash / pu_;
        down_ = (ermqdt - pu_ * up_) / (1.0 - pu_);
    }

    Real ExtendedJoshi43::underlying(Size i, Size index) const {
        Time stepTime = i*this->dt_;
        Real variance = this->treeProcess_->variance(stepTime, x0_, end_);

        Real driftStep_ = this->driftStepStorage[i];

        Real ermqdt = std::exp(driftStep_ + 0.5*variance / oddSteps_);
		Real d2 = (std::log(x0_ / strike_) + driftStep_*oddSteps_) /
			std::sqrt(variance);

        Real pu = computeUpProb((oddSteps_-1.0)/2.0,d2 );
        Real pdash = computeUpProb((oddSteps_-1.0)/2.0,d2+std::sqrt(variance));
        Real up = ermqdt * pdash / pu;
        Real down = (ermqdt - pu * up) / (1.0 - pu);

        return x0_ * std::pow(down, Real(BigInteger(i)-BigInteger(index)))
            * std::pow(up, Real(index));
    }

    Real ExtendedJoshi43::probability(Size i, Size, Size branch) const {
        Time stepTime = i*this->dt_;
        Real variance = this->treeProcess_->variance(stepTime, x0_, end_);

        Real d2 = (std::log(x0_ / strike_) + this->driftStepStorage[i]*oddSteps_) /
			std::sqrt(variance);

        Real pu = computeUpProb((oddSteps_-1.0)/2.0,d2 );
        Real pd = 1.0 - pu;

        return (branch == 1 ? pu : pd);
    }

}

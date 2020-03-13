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

/*! \file extendedbinomialtree.hpp
    \brief Time-dependent binomial tree class
*/

#ifndef quantlib_extended_binomial_tree3_hpp
#define quantlib_extended_binomial_tree3_hpp
//#include "extendedbinomialtree3.cpp"
#include <ql/instruments/dividendschedule.hpp>
#include <ql/methods/lattices/tree.hpp>
#include <ql/stochasticprocess.hpp>
#include <vector>
#include <ql/functional.hpp>

namespace QuantLib {

    //! Binomial tree base class
    /*! \ingroup lattices */
    //using namespace ext::placeholders;

    template <class T>
    class ExtendedBinomialTree3 : public Tree<T> {
      public:
        enum Branches { branches = 2 };
        ExtendedBinomialTree3(
                        const ext::shared_ptr<StochasticProcess1D>& process,
                        Time end,
                        Size steps)
        : Tree<T>(steps+1), treeProcess_(process) {
            x0_ = process->x0();
            dt_ = end/steps;
            driftPerStep_ = process->drift(0.0, x0_) * dt_;
            for (Size i = 0; i <= steps; i ++) {
                driftStepStorage.push_back(this->driftStep(i));
            }

        }
        Size size(Size i) const {
            return i+1;
        }
        Size descendant(Size, Size index, Size branch) const {
            return index + branch;
        }
      protected:
        Real driftStep(Time driftTime) const {
            return this->treeProcess_->drift(driftTime, x0_) * dt_;
        }
        std::vector<Real> driftStepStorage;
        Real x0_, driftPerStep_;
        Time dt_;

      protected:
        ext::shared_ptr<StochasticProcess1D> treeProcess_;
    };


    //! Base class for equal probabilities binomial tree
    /*! \ingroup lattices */
    template <class T>
    class ExtendedEqualProbabilitiesBinomialTree3
        : public ExtendedBinomialTree3<T> {
      public:
        ExtendedEqualProbabilitiesBinomialTree3(
                        const ext::shared_ptr<StochasticProcess1D>& process,
                        Time end,
                        Size steps)
        : ExtendedBinomialTree3<T>(process, end, steps) {}
        virtual ~ExtendedEqualProbabilitiesBinomialTree3() {}

        Real underlying(Size i, Size index) const {
            BigInteger j = 2*BigInteger(index) - BigInteger(i);
            return this->x0_*std::exp(i*this->driftStepStorage[i] + j*this->upStepStorage[i]);
        }

        Real probability(Size, Size, Size) const { return 0.5; }
      protected:
        virtual Real upStep(Size i) const = 0;
        std::vector<Real> upStepStorage;
        Real up_;
    };


    //! Base class for equal jumps binomial tree
    /*! \ingroup lattices */
    template <class T>
    class ExtendedEqualJumpsBinomialTree3 : public ExtendedBinomialTree3<T> {
      public:
        ExtendedEqualJumpsBinomialTree3(
                        const ext::shared_ptr<StochasticProcess1D>& process,
                        Time end,
                        Size steps)
        : ExtendedBinomialTree3<T>(process, end, steps) {}
        virtual ~ExtendedEqualJumpsBinomialTree3() {}
        Real underlying(Size i, Size index) const {
            BigInteger j = 2*BigInteger(index) - BigInteger(i);
            return this->x0_*std::exp(j*this->dxStepStorage[i]);

        }

        Real probability(Size i, Size, Size branch) const {
            Real upProb = this->probUpStorage[i];
            Real downProb = 1 - upProb;
            return (branch == 1 ? upProb : downProb);
        }
      protected:
        //probability of a up move
        virtual Real probUp(Size i) const = 0;
        //time dependent term dx_
        virtual Real dxStep(Size i) const = 0;
        std::vector<Real> probUpStorage;
        std::vector<Real> dxStepStorage;
        Real dx_, pu_, pd_;
    };


    //! Jarrow-Rudd (multiplicative) equal probabilities binomial tree
    /*! \ingroup lattices */
    class ExtendedJarrowRudd3
        : public ExtendedEqualProbabilitiesBinomialTree3<ExtendedJarrowRudd3> {
      public:
        ExtendedJarrowRudd3(const ext::shared_ptr<StochasticProcess1D>&,
                           Time end,
                           Size steps,
                           Real strike);
      protected:
        Real upStep(Size i) const;
    };


    //! Cox-Ross-Rubinstein (multiplicative) equal jumps binomial tree
    /*! \ingroup lattices */
    class ExtendedCoxRossRubinstein3
        : public ExtendedEqualJumpsBinomialTree3<ExtendedCoxRossRubinstein3> {
      public:
        ExtendedCoxRossRubinstein3(const ext::shared_ptr<StochasticProcess1D>&,
                                  Time end,
                                  Size steps,
                                  Real strike);
      protected:
          Real probUp(Size i) const;
          Real dxStep(Size i) const;
    };


    //! Additive equal probabilities binomial tree
    /*! \ingroup lattices */
    class ExtendedAdditiveEQPBinomialTree3
        : public ExtendedEqualProbabilitiesBinomialTree3<
                                            ExtendedAdditiveEQPBinomialTree3> {
      public:
        ExtendedAdditiveEQPBinomialTree3(
                        const ext::shared_ptr<StochasticProcess1D>&,
                        Time end,
                        Size steps,
                        Real strike);

      protected:
          Real upStep(Size i) const;
    };


    //! %Trigeorgis (additive equal jumps) binomial tree
    /*! \ingroup lattices */
    class ExtendedTrigeorgis3
        : public ExtendedEqualJumpsBinomialTree3<ExtendedTrigeorgis3> {
      public:
        ExtendedTrigeorgis3(const ext::shared_ptr<StochasticProcess1D>&,
                           Time end,
                           Size steps,
                           Real strike);
    protected:
        Real probUp(Size i) const;
        Real dxStep(Size i) const;
    };


    //! %Tian tree: third moment matching, multiplicative approach
    /*! \ingroup lattices */
    class ExtendedTian3 : public ExtendedBinomialTree3<ExtendedTian3> {
      public:
        ExtendedTian3(const ext::shared_ptr<StochasticProcess1D>&,
                     Time end,
                     Size steps,
                     Real strike);

        Real underlying(Size i, Size index) const;
        Real probability(Size, Size, Size branch) const;
      protected:
        Real up_, down_, pu_, pd_;
    };

    //! Leisen & Reimer tree: multiplicative approach
    /*! \ingroup lattices */
    class ExtendedLeisenReimer3
        : public ExtendedBinomialTree3<ExtendedLeisenReimer3> {
      public:
        ExtendedLeisenReimer3(const ext::shared_ptr<StochasticProcess1D>&,
                             Time end,
                             Size steps,
                             Real strike);

        Real underlying(Size i, Size index) const;
        Real probability(Size, Size, Size branch) const;
      protected:
        Time end_;
        Size oddSteps_;
        Real strike_, up_, down_, pu_, pd_;
    };


     class ExtendedJoshi43 : public ExtendedBinomialTree3<ExtendedJoshi43> {
      public:
        ExtendedJoshi43(const ext::shared_ptr<StochasticProcess1D>&,
                       Time end,
                       Size steps,
                       Real strike);

        Real underlying(Size i, Size index) const;
        Real probability(Size, Size, Size branch) const;
      protected:
        Real computeUpProb(Real k, Real dj) const;
        Time end_;
        Size oddSteps_;
        Real strike_, up_, down_, pu_, pd_;
    };


}


#endif

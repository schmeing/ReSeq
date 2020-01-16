#ifndef REFERENCETEST_H
#define REFERENCETEST_H
#include "Reference.h"

#include <stdint.h>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class ReferenceTest : public BasicTestClass{
	public:
		static void Register();

	protected:
		Reference ref_;

		void TestVariantClass();
		void TestInsertVariant();
		void TestVariationLoading();
		void TestVariationPositionLoading();
		void TestLoadingAndAccess();
		void TestGC();
		void TestSumBias();
		void TestGetFragmentSites();
		void TestReplaceN();
		void TestExclusionRegions();
		void TestMethylationLoading();
	};
}

#endif // REFERENCETEST_H

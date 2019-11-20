#ifndef ADAPTERSTATSTEST_H
#define ADAPTERSTATSTEST_H
#include "AdapterStats.h"

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class AdapterStatsTest : public BasicTestClass{
	public:
		static void Register();

	protected:
		AdapterStats *test_;

		void CreateTestObject();
		void DeleteTestObject();

		void LoadAdapters();
		static void TestLoading(const AdapterStats &test);
		static void TestSumming(AdapterStats &test);

		virtual void TearDown();
	public:
		AdapterStatsTest():
			test_(NULL)
			{
		}

		static void TestAdapters(const AdapterStats &test, const char *context);
		static void TestNexteraAdapters(const AdapterStats &test);
	};
}

#endif // ADAPTERSTATSTEST_H

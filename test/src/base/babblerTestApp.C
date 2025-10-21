//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "babblerTestApp.h"
#include "babblerApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
babblerTestApp::validParams()
{
  InputParameters params = babblerApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

babblerTestApp::babblerTestApp(const InputParameters & parameters) : MooseApp(parameters)
{
  babblerTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

babblerTestApp::~babblerTestApp() {}

void
babblerTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  babblerApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"babblerTestApp"});
    Registry::registerActionsTo(af, {"babblerTestApp"});
  }
}

void
babblerTestApp::registerApps()
{
  registerApp(babblerApp);
  registerApp(babblerTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
babblerTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  babblerTestApp::registerAll(f, af, s);
}
extern "C" void
babblerTestApp__registerApps()
{
  babblerTestApp::registerApps();
}

// -*-c++-*-
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#include <masa_internal.h>
#include <cmath>

using namespace MASA;
using namespace std;

template <typename Scalar>
MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::navierstokes_3d_incompressible_homogeneous()
{
  this->mmsname = "navierstokes_3d_incompressible_homogeneous";
  this->dimension = 3;

  this->register_var("a",&a);
  this->register_var("b",&b);
  this->register_var("c",&c);
  this->register_var("d",&d);
  this->register_var("beta",&beta);
  this->register_var("gamma",&gamma);
  this->register_var("delta",&delta);
  this->register_var("nu",&nu);
  this->register_var("kx",&kx);
  this->register_var("kz",&kz);
  this->register_var("ky",&ky);

  this->init_var();

} // done with constructor

template <typename Scalar>
int MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("a",1.05);
  err += this->set_var("b",2.15);
  err += this->set_var("c",-3.2);
  err += this->set_var("d",10.1);
  err += this->set_var("beta",2.2);
  err += this->set_var("gamma",2.4);
  err += this->set_var("delta",2.0);
  err += this->set_var("nu",.02);
  err += this->set_var("kx",1);
  err += this->set_var("kz",1);
  err += this->set_var("ky",1);

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------

// u component of velocity source term
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::eval_q_u(Scalar x, Scalar y, Scalar z)
{
  // NS equation residuals
  Scalar Q_u = 
    (2.*pow(a,2)*kx*pow(ky,2)*pow(kz,2)*cos(kx*x)*pow(cos(ky*y),2)*
      pow(cos(kz*z),2))/
    (pow(beta + sin(kx*x),3)*pow(delta + sin(ky*y),4)*
      pow(gamma + sin(kz*z),4)) + 
   (3.*a*b*kx*pow(ky,2)*pow(kz,2)*cos(kx*x)*pow(cos(ky*y),2)*
      pow(cos(kz*z),2))/
    (pow(beta + sin(kx*x),3)*pow(delta + sin(ky*y),4)*
      pow(gamma + sin(kz*z),4)) + 
   (3.*a*c*kx*pow(ky,2)*pow(kz,2)*cos(kx*x)*pow(cos(ky*y),2)*
      pow(cos(kz*z),2))/
    (pow(beta + sin(kx*x),3)*pow(delta + sin(ky*y),4)*
      pow(gamma + sin(kz*z),4)) + 
   (1.*a*b*kx*pow(ky,2)*pow(kz,2)*cos(kx*x)*pow(cos(kz*z),2)*
      sin(ky*y))/
    (pow(beta + sin(kx*x),3)*pow(delta + sin(ky*y),3)*
      pow(gamma + sin(kz*z),4)) + 
   (1.*a*c*kx*pow(ky,2)*pow(kz,2)*cos(kx*x)*pow(cos(ky*y),2)*
      sin(kz*z))/
    (pow(beta + sin(kx*x),3)*pow(delta + sin(ky*y),4)*
      pow(gamma + sin(kz*z),3)) + 
   (1.*d*kx*cos(kx*x))/
    (pow(beta + sin(kx*x),2)*(delta + sin(ky*y))*(gamma + sin(kz*z))) + 
   nu*((6.*a*ky*pow(kz,3)*cos(ky*y)*pow(cos(kz*z),3))/
       ((beta + sin(kx*x))*pow(delta + sin(ky*y),2)*
         pow(gamma + sin(kz*z),4)) + 
      (6.*a*ky*pow(kz,3)*cos(ky*y)*cos(kz*z)*sin(kz*z))/
       ((beta + sin(kx*x))*pow(delta + sin(ky*y),2)*
         pow(gamma + sin(kz*z),3)) + 
      (6.*a*pow(ky,3)*kz*pow(cos(ky*y),3)*cos(kz*z))/
       ((beta + sin(kx*x))*pow(delta + sin(ky*y),4)*
         pow(gamma + sin(kz*z),2)) + 
      (6.*a*pow(ky,3)*kz*cos(ky*y)*cos(kz*z)*sin(ky*y))/
       ((beta + sin(kx*x))*pow(delta + sin(ky*y),3)*
         pow(gamma + sin(kz*z),2)) + 
      (2.*a*pow(kx,2)*ky*kz*pow(cos(kx*x),2)*cos(ky*y)*cos(kz*z))/
       (pow(beta + sin(kx*x),3)*pow(delta + sin(ky*y),2)*
         pow(gamma + sin(kz*z),2)) + 
      (1.*a*pow(kx,2)*ky*kz*cos(ky*y)*cos(kz*z)*sin(kx*x))/
       (pow(beta + sin(kx*x),2)*pow(delta + sin(ky*y),2)*
         pow(gamma + sin(kz*z),2)) - 
      (1.*a*pow(ky,3)*kz*cos(ky*y)*cos(kz*z))/
       ((beta + sin(kx*x))*pow(delta + sin(ky*y),2)*
         pow(gamma + sin(kz*z),2)) - 
      (1.*a*ky*pow(kz,3)*cos(ky*y)*cos(kz*z))/
       ((beta + sin(kx*x))*pow(delta + sin(ky*y),2)*
         pow(gamma + sin(kz*z),2)));

  return -Q_u;
}

// v component of velocity source term
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::eval_q_v(Scalar x, Scalar y, Scalar z)
{
    Scalar Q_v = 
    (3.*a*b*pow(kx,2)*ky*pow(kz,2)*pow(cos(kx*x),2)*cos(ky*y)*
      pow(cos(kz*z),2))/
    (pow(beta + sin(kx*x),4)*pow(delta + sin(ky*y),3)*
      pow(gamma + sin(kz*z),4)) + 
   (2.*pow(b,2)*pow(kx,2)*ky*pow(kz,2)*pow(cos(kx*x),2)*cos(ky*y)*
      pow(cos(kz*z),2))/
    (pow(beta + sin(kx*x),4)*pow(delta + sin(ky*y),3)*
      pow(gamma + sin(kz*z),4)) + 
   (3.*b*c*pow(kx,2)*ky*pow(kz,2)*pow(cos(kx*x),2)*cos(ky*y)*
      pow(cos(kz*z),2))/
    (pow(beta + sin(kx*x),4)*pow(delta + sin(ky*y),3)*
      pow(gamma + sin(kz*z),4)) + 
   (1.*a*b*pow(kx,2)*ky*pow(kz,2)*cos(ky*y)*pow(cos(kz*z),2)*
      sin(kx*x))/
    (pow(beta + sin(kx*x),3)*pow(delta + sin(ky*y),3)*
      pow(gamma + sin(kz*z),4)) + 
   (1.*b*c*pow(kx,2)*ky*pow(kz,2)*pow(cos(kx*x),2)*cos(ky*y)*
      sin(kz*z))/
    (pow(beta + sin(kx*x),4)*pow(delta + sin(ky*y),3)*
      pow(gamma + sin(kz*z),3)) + 
   (1.*d*ky*cos(ky*y))/
    ((beta + sin(kx*x))*pow(delta + sin(ky*y),2)*(gamma + sin(kz*z))) + 
   nu*((6.*b*kx*pow(kz,3)*cos(kx*x)*pow(cos(kz*z),3))/
       (pow(beta + sin(kx*x),2)*(delta + sin(ky*y))*
         pow(gamma + sin(kz*z),4)) + 
      (6.*b*kx*pow(kz,3)*cos(kx*x)*cos(kz*z)*sin(kz*z))/
       (pow(beta + sin(kx*x),2)*(delta + sin(ky*y))*
         pow(gamma + sin(kz*z),3)) + 
      (2.*b*kx*pow(ky,2)*kz*cos(kx*x)*pow(cos(ky*y),2)*cos(kz*z))/
       (pow(beta + sin(kx*x),2)*pow(delta + sin(ky*y),3)*
         pow(gamma + sin(kz*z),2)) + 
      (1.*b*kx*pow(ky,2)*kz*cos(kx*x)*cos(kz*z)*sin(ky*y))/
       (pow(beta + sin(kx*x),2)*pow(delta + sin(ky*y),2)*
         pow(gamma + sin(kz*z),2)) + 
      (6.*b*pow(kx,3)*kz*pow(cos(kx*x),3)*cos(kz*z))/
       (pow(beta + sin(kx*x),4)*(delta + sin(ky*y))*
         pow(gamma + sin(kz*z),2)) + 
      (6.*b*pow(kx,3)*kz*cos(kx*x)*cos(kz*z)*sin(kx*x))/
       (pow(beta + sin(kx*x),3)*(delta + sin(ky*y))*
         pow(gamma + sin(kz*z),2)) - 
      (1.*b*pow(kx,3)*kz*cos(kx*x)*cos(kz*z))/
       (pow(beta + sin(kx*x),2)*(delta + sin(ky*y))*
         pow(gamma + sin(kz*z),2)) - 
      (1.*b*kx*pow(kz,3)*cos(kx*x)*cos(kz*z))/
       (pow(beta + sin(kx*x),2)*(delta + sin(ky*y))*
         pow(gamma + sin(kz*z),2)));

  return -Q_v;

}

// w component of velocity source term
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::eval_q_w(Scalar x, Scalar y, Scalar z)
{
  Scalar Q_w = 
    (3.*a*c*pow(kx,2)*pow(ky,2)*kz*pow(cos(kx*x),2)*pow(cos(ky*y),2)*
      cos(kz*z))/
    (pow(beta + sin(kx*x),4)*pow(delta + sin(ky*y),4)*
      pow(gamma + sin(kz*z),3)) + 
   (3.*b*c*pow(kx,2)*pow(ky,2)*kz*pow(cos(kx*x),2)*
      pow(cos(ky*y),2)*cos(kz*z))/
    (pow(beta + sin(kx*x),4)*pow(delta + sin(ky*y),4)*
      pow(gamma + sin(kz*z),3)) + 
   (2.*pow(c,2)*pow(kx,2)*pow(ky,2)*kz*pow(cos(kx*x),2)*
      pow(cos(ky*y),2)*cos(kz*z))/
    (pow(beta + sin(kx*x),4)*pow(delta + sin(ky*y),4)*
      pow(gamma + sin(kz*z),3)) + 
   (1.*a*c*pow(kx,2)*pow(ky,2)*kz*pow(cos(ky*y),2)*cos(kz*z)*
      sin(kx*x))/
    (pow(beta + sin(kx*x),3)*pow(delta + sin(ky*y),4)*
      pow(gamma + sin(kz*z),3)) + 
   (1.*b*c*pow(kx,2)*pow(ky,2)*kz*pow(cos(kx*x),2)*cos(kz*z)*
      sin(ky*y))/
    (pow(beta + sin(kx*x),4)*pow(delta + sin(ky*y),3)*
      pow(gamma + sin(kz*z),3)) + 
   (1.*d*kz*cos(kz*z))/
    ((beta + sin(kx*x))*(delta + sin(ky*y))*pow(gamma + sin(kz*z),2)) + 
   nu*((2.*c*kx*ky*pow(kz,2)*cos(kx*x)*cos(ky*y)*pow(cos(kz*z),2))/
       (pow(beta + sin(kx*x),2)*pow(delta + sin(ky*y),2)*
         pow(gamma + sin(kz*z),3)) + 
      (1.*c*kx*ky*pow(kz,2)*cos(kx*x)*cos(ky*y)*sin(kz*z))/
       (pow(beta + sin(kx*x),2)*pow(delta + sin(ky*y),2)*
         pow(gamma + sin(kz*z),2)) + 
      (6.*c*kx*pow(ky,3)*cos(kx*x)*pow(cos(ky*y),3))/
       (pow(beta + sin(kx*x),2)*pow(delta + sin(ky*y),4)*
         (gamma + sin(kz*z))) + 
      (6.*c*kx*pow(ky,3)*cos(kx*x)*cos(ky*y)*sin(ky*y))/
       (pow(beta + sin(kx*x),2)*pow(delta + sin(ky*y),3)*
         (gamma + sin(kz*z))) + 
      (6.*c*pow(kx,3)*ky*pow(cos(kx*x),3)*cos(ky*y))/
       (pow(beta + sin(kx*x),4)*pow(delta + sin(ky*y),2)*
         (gamma + sin(kz*z))) + 
      (6.*c*pow(kx,3)*ky*cos(kx*x)*cos(ky*y)*sin(kx*x))/
       (pow(beta + sin(kx*x),3)*pow(delta + sin(ky*y),2)*
         (gamma + sin(kz*z))) - 
      (1.*c*pow(kx,3)*ky*cos(kx*x)*cos(ky*y))/
       (pow(beta + sin(kx*x),2)*pow(delta + sin(ky*y),2)*
         (gamma + sin(kz*z))) - 
      (1.*c*kx*pow(ky,3)*cos(kx*x)*cos(ky*y))/
       (pow(beta + sin(kx*x),2)*pow(delta + sin(ky*y),2)*
         (gamma + sin(kz*z))));
         
  return -Q_w;
  
}



// ----------------------------------------
// Analytical Terms
// ----------------------------------------

//
// main functions
// 

// example of a public method called from eval_exact_t
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::eval_exact_u(Scalar x, Scalar y, Scalar z)
{
  Scalar exact_u =   (a*ky*kz*cos(ky*y)*cos(kz*z))/
   ((beta + sin(kx*x))*pow(delta + sin(ky*y),2)*pow(gamma + sin(kz*z),2));
  return exact_u;
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::eval_exact_v(Scalar x, Scalar y, Scalar z)
{
  Scalar exact_v = (b*kx*kz*cos(kx*x)*cos(kz*z))/
   (pow(beta + sin(kx*x),2)*(delta + sin(ky*y))*pow(gamma + sin(kz*z),2));
  return exact_v;
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::eval_exact_w(Scalar x, Scalar y, Scalar z)
{
  Scalar exact_w = (c*kx*ky*cos(kx*x)*cos(ky*y))/
   (pow(beta + sin(kx*x),2)*pow(delta + sin(ky*y),2)*(gamma + sin(kz*z)));
  return exact_w;
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_incompressible_homogeneous<Scalar>::eval_exact_p(Scalar x, Scalar y, Scalar z)
{
  Scalar P = (d)/((beta + sin(kx*x))*(delta + sin(ky*y))*(gamma + sin(kz*z)));
  return P;
}


// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::navierstokes_3d_incompressible_homogeneous);



//---------------------------------------------------------
// AUTOMASA
// Generated on: 2013-05-08 11:32:28
//---------------------------------------------------------

//void FVHMComponent::calculateCellVelocityGradients()
//{
//#ifdef USE_OPENMP
//#pragma omp parallel for
//#endif
//  for(int i = 0 ; i < m_controlVolumes.size() ; i++)
//  {
//    TriCV *cv = m_controlVolumes[i];

//    if(cv->isActive)
//    {
//      //too costly
//      Vect v1 = cv->e_n[0] * cv->edgeNormVels[0].Value;
//      Vect v2 = cv->e_n[1] * cv->edgeNormVels[1].Value;
//      Vect v3 = cv->e_n[2] * cv->edgeNormVels[2].Value;

//      double values[3];
//      values[0] = v1.x;
//      values[1] = v2.x;
//      values[2] = v3.x;

//      Vect grad = TriCV::lsGradReconstruction(cv->vels[0].Value , cv->r_e, values , 3);

//      cv->grad_u->x = grad.x;
//      cv->grad_u->y = grad.y;

//      values[0] = v1.y;
//      values[1] = v2.y;
//      values[2] = v3.y;

//      grad = TriCV::lsGradReconstruction(cv->vels[1].Value , cv->r_e, values, 3);

//      cv->grad_v->x = grad.x ;
//      cv->grad_v->y = grad.y ;

//    }
//  }
//}
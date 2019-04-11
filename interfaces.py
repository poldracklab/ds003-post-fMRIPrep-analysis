#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Additional interfaces required for workflow. To be upstreamed into Nipype.
"""
from nipype.interfaces.fsl.base import FSLCommand, FSLCommandInputSpec
from nipype.interfaces.base import traits, TraitedSpec


class PtoZInputSpec(FSLCommandInputSpec):
    pvalue = traits.Float(0.05, argstr='%f', mandatory=True, usedefault=True,
                          position=0,
                          desc='p-value for which the corresponding '
                               'z-statistic should be computed')
    two_tailed = traits.Bool(argstr='-2', position=1,
                             desc='use 2-tailed conversion (default is '
                                  '1-tailed)')
    grf = traits.Float(argstr='-g %f', position=2,
                       desc='use GRF maximum-height theory instead of '
                            'Gaussian PDF. To enable this option, specify '
                            'the number of resels as the argument. This can '
                            'be estimated using fsl.SmoothEstimate.')


class PtoZOutputSpec(TraitedSpec):
    zstat = traits.Float(
        desc='z-statistic corresponding to specified p-value')


class PtoZ(FSLCommand):
    """Determine the z-value corresponding to an observed p-value."""
    input_spec = PtoZInputSpec
    output_spec = PtoZOutputSpec
    _cmd = 'ptoz'

    def aggregate_outputs(self, runtime=None, needed_outputs=None):
        outputs = self._outputs()
        outputs.zstat = float(runtime.stdout)
        return outputs

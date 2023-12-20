.. currentmodule:: pywatershed

####################################################
pywatershed: A hydrologic model in Python
####################################################

Welcome to the `pywatershed` docs!

Pywatershed is Python package for simulating hydrologic processes motivated by
the need to modernize important, legacy hydrologic models at the USGS,
particularly the
`Precipitation-Runoff Modeling System <https://www.usgs.gov/software/precipitation-runoff-modeling-system-prms>`__
(PRMS, Markstrom et al., 2015) and its role in
`GSFLOW <https://www.usgs.gov/software/gsflow-coupled-groundwater-and-surface-water-flow-model>`__
(Markstrom et al., 2008).
The goal of modernization is to make these legacy models more flexible as process
representations, to support testing of alternative hydrologic process
conceptualizations, and to facilitate the incorporation of cutting edge
modeling techniques and data sources. Pywatershed is a place for experimentation
with software design, process representation, and data fusion in the context
of well-established hydrologic process modeling.

The Python language was choosen because it is accessible to a wide audience of
potential contributors which will help foster community development and
experimentation. A large number of advanced libraries available for Python can
also be applied to hdyrologic modeling, including libraries for parallelism, data
access and manipulation, and machine learning.

Following the conceptual design of PRMS, pywatershed calculates explicit solutions
of spatially distributed hydrologic process representations including evaporation,
transpiration, runoff, infiltration, interflow, snowpack, soil moisture, conceptual
groundwater storage, and channel flow. These process representations simulate hydrologic
response and water budgets given inputs of spatially distributed weather variables and
land use change at temporal scales ranging from days to centuries. 

Pywatershed enhances PRMS with a new software design that is object-oriented and highly
flexible, allowing users to easily run "sub-models", replace process representations, and
incorporate new data. There are base classes which manage mass and energy conservation
and the implementation of concrete process classes follows a self-describing design
which allows for a Model class to connect hydrologic process classes based on their
descriptions of themselves. A variety of input data sources is managed by the
Adapter class which implements subclasses for different sources. The design of
pywatershed is documented in these docs and also demonstrated by numbered Jupyter
Notebooks in the `examples/` directory.

The flexible structure of pywatershed helps it to couple with other hydrologic
models. We can easily one-way couple pywatershed to
`MODFLOW 6 <https://www.usgs.gov/software/modflow-6-usgs-modular-hydrologic-model>`__
(MF6, Hughes et al., 2017)
via its
`XMI interface <https://www.usgs.gov/publications/modflow-application-programming-interface-simulationcontrol-and-software>`__
(Hughes et al., 2022).
We are working towards a two-way, tight coupling with MF6 to reproduce GSFLOW. Our goal is
support integrated hydrologic process modeling of surface water and groundwater in a
sustainable manner that allows individual software components to evolve independently. 


=========================
Current version: 1.0.0
=========================
With pywatershed version 1.0.0, we have faithfully reproduced the PRMS process representations used in
the USGS `National Hydrolgical Model <https://pubs.usgs.gov/publication/tm6B9>`__ (NHM, Regan et al.,
2018). For more information on version 1.0.0 see the
`release notes <https://github.com/EC-USGS/pywatershed/releases/tag/1.0.0>`_
and the `extended release notes <https://ec-usgs.github.io/pywatershed/2023/12/18/v1-0-0-overview>`_
for version 1.0.0.

============================================
Upcoming development in 2024
============================================
The broad goal is to reproduce GSFLOW coupling using the MODFLOW 6 API. This will include
gridded configurations and cascading flows.
We are also working on reservoir representations.

=================
Getting started
=================
Please note that you can browse the API reference, developer info, and index
using the table of contents on the left. But *the best way to get started
with pywatershed is to dive into the example notebooks*.

| For introductory example notebooks, look in the `examples/ <https://github.com/EC-USGS/pywatershed/tree/main/examples>`_ directory in the repository. Numbered starting at 00, these are meant to be completed in order. Notebook outputs are not saved in Github. But you can run these notebooks locally or using WholeTale (free but sign-up or log-in required) where the pywatershed environment is all ready to go:

.. image:: https://raw.githubusercontent.com/whole-tale/wt-design-docs/master/badges/wholetale-explore.svg
   :target: https://dashboard.wholetale.org

* `Run the latest release in WholeTale <https://dashboard.wholetale.org/run/64ae29e8a887f48b9f173678?tab=metadata>`_
* `Run the develop branch in WholeTale <https://dashboard.wholetale.org/run/64ae25c3a887f48b9f1735c8?tab=metadata>`_

See `README.md <https://github.com/EC-USGS/pywatershed/tree/develop/README.md>`_ for more details
on both `running locally <https://github.com/EC-USGS/pywatershed#installation>`_
or `using WholeTale <https://github.com/EC-USGS/pywatershed#example-notebooks>`_.

========================
Community engagement
========================
We value your feedback! Please use `discussions <https://github.com/EC-USGS/pywatershed/discussions>`_
or `issues <https://github.com/EC-USGS/pywatershed/issues>`_ on Github. You may also suggest
edits to this documentation or open an issue by clicking on the Github Octocat at the top of the page.
For more in-depth contributions, please start by reading over
the `DEVELOPER.md file <https://github.com/EC-USGS/pywatershed/blob/develop/DEVELOPER.md>`_.

Thank you for your interest.

===========
References
===========
* `Hughes, J. D., Langevin, C. D., & Banta, E. R. (2017). Documentation for the MODFLOW 6 framework (No. 6-A57). US Geological Survey. <https://pubs.usgs.gov/publication/tm6A57>`__
* `Hughes, J. D., Russcher, M. J., Langevin, C. D., Morway, E. D., & McDonald, R. R. (2022). The MODFLOW Application Programming Interface for simulation control and software interoperability. Environmental Modelling & Software, 148, 105257. <https://www.sciencedirect.com/science/article/pii/S1364815221002991>`__
* `Markstrom, S. L., Niswonger, R. G., Regan, R. S., Prudic, D. E., & Barlow, P. M. (2008). GSFLOW-Coupled Ground-water and Surface-water FLOW model based on the integration of the Precipitation-Runoff Modeling System (PRMS) and the Modular Ground-Water Flow Model (MODFLOW-2005). US Geological Survey techniques and methods, 6, 240. <https://pubs.usgs.gov/tm/tm6d1/>`__
* `Markstrom, S. L., Regan, R. S., Hay, L. E., Viger, R. J., Webb, R. M., Payn, R. A., & LaFontaine, J. H. (2015). PRMS-IV, the precipitation-runoff modeling system, version 4 (No. 6-B7). US Geological Survey. <https://pubs.usgs.gov/tm/6b7/>`__
* `Regan, R. S., Markstrom, S. L., Hay, L. E., Viger, R. J., Norton, P. A., Driscoll, J. M., & LaFontaine, J. H. (2018). Description of the national hydrologic model for use with the precipitation-runoff modeling system (prms) (No. 6-B9). US Geological Survey. <https://pubs.usgs.gov/publication/tm6B9>`__


.. toctree::
   :caption: Home

   self
	     

  
.. toctree::
   :hidden:
   :caption: API Reference

    API Summary <api.rst>
    meta <api/generated/pywatershed.meta.rst>
    Control <api/generated/pywatershed.Control.rst>
    Parameters <api/parameters>
    adapter <api/adapter>
    Atmosphere <api/atmosphere>
    Hydrology <api/hydrology>
    Model <api/generated/pywatershed.Model.rst>
    Base Classes <api/base>
    Utilities <api/utils>

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: For developers

   What’s New <whats-new>
   GitHub repository <https://github.com/EC-USGS/pywatershed>

.. toctree::
   :hidden:
   :caption: Index

    Index <genindex>
..
    Module Index <modindex>

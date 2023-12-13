.. currentmodule:: pywatershed

####################################################
pywatershed: A hydrologic model in Python
####################################################

Welcome to the `pywatershed` docs!

Pywatershed is a path for modernizing legacy USGS hydrology models,
in particular
`PRMS <https://www.usgs.gov/software/precipitation-runoff-modeling-system-prms>`__
and its role in
`GSFLOW <https://www.usgs.gov/software/gsflow-coupled-groundwater-and-surface-water-flow-model>`__.
Pywatershed is also a place for experimentation with these hydrologic models in this context.


Following PRMS, pywatershed calculates explicit solutions of spatially
distributed hydrologic process representations including evaporation, transpiration,
runoff, infiltration, interflow, snowpack, soil moisture, conceptual groundwater storage,
and channel flow. These process representations simulate hydrologic response and water
budgets given inputs of spatially distributed weather variables and land use change at
temporal scales ranging from days to centuries. 

Pywatershed enhances PRMS with a new software design that is object-oriented and highly
flexible, allowing users to easily run "sub-models", replace process representations, and
incorporate new data. The Python language is more accessible to a wider audience of
potential contributors which can help foster community development. A large number of
advanced libraries are available within Python that can be used for experimentation,
including libraries for parallelism, data access and manipulation, and machine learning.
Finally, we can easily couple pywatershed to
`MODFLOW 6 <https://www.usgs.gov/software/modflow-6-usgs-modular-hydrologic-model>`__
via its
`XMI interface <https://www.usgs.gov/publications/modflow-application-programming-interface-simulationcontrol-and-software>`__.
This coupling will eventually bring pywatershed and MODFLOW 6 together to reproduce
GSFLOW using a modern software design.


=========================
Current Version: 1.0.0
=========================
With pywatershed version 1.0.0, we have faithfully reproduced the NHM
configuration of PRMS. For more information see the
`release notes <https://github.com/EC-USGS/pywatershed/releases/tag/1.0.0>`_
and the `extended release notes <https://ec-usgs.github.io/pywatershed/2023/11/14/v1-0-0-overview>`_
for version 1.0.0.

============================================
Upcoming development in 2024
============================================
The broad goal is to reproduce GSFLOW coupling using the MODFLOW 6 API. This will include
gridded configurations and cascading flows.
We are also working on reservoir representations.

=================
Getting Started
=================
Please note that you can browse the API reference, developer info, and index
using the table of contents on the left. But *the best way to get started
with pywatershed is to dive in to the example notebooks*.

| For introductory example notebooks, look in the `examples/ <https://github.com/EC-USGS/pywatershed/tree/main/examples>`_ directory in the repository. Numbered starting at 00, these are meant to be completed in order. Notebook outputs are not saved in Github. But you can run these notebooks locally or using WholeTale (free but sign-up or log-in required) where the pywatershed environment is all ready to go:

.. image:: https://raw.githubusercontent.com/whole-tale/wt-design-docs/master/badges/wholetale-explore.svg
   :target: https://dashboard.wholetale.org

* `Run the latest release in WholeTale <https://dashboard.wholetale.org/run/64ae29e8a887f48b9f173678?tab=metadata>`_
* `Run the develop branch in WholeTale <https://dashboard.wholetale.org/run/64ae25c3a887f48b9f1735c8?tab=metadata>`_

See `README.md <https://github.com/EC-USGS/pywatershed/tree/develop/README.md>`_ for more details
on both `running locally <https://github.com/EC-USGS/pywatershed#installation>`_
or `using WholeTale <https://github.com/EC-USGS/pywatershed#example-notebooks>`_.

========================
Community Engagement
========================
We value your feedback! Please use `discussions <https://github.com/EC-USGS/pywatershed/discussions>`_
or `issues <https://github.com/EC-USGS/pywatershed/issues>`_ on our Github page. You may also suggest
edits to this documentation or open an issue by clicking on the Github Octocat at the top of the page.
For more in-depth contributions, please start by reading over
the `DEVELOPER.md file <https://github.com/EC-USGS/pywatershed/blob/develop/DEVELOPER.md>`_.

Thank you for your interest.


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

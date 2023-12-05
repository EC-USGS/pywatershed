.. currentmodule:: pywatershed

####################################################
pywatershed: a hydrologic model in Python
####################################################

Welcome to the `pywatershed` docs!

pywatershed is a path for modernizing legacy USGS hydrology models,
in particular PRMS and GSFlow. It is also a place for experimentation with
hydrologic models in this context.

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

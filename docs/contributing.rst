=======================
Contributing Guidelines
=======================

ChemGED was a side project that came about to support separate research on semantically aware
representation learning for chemical graphs. As a result, we only ever implemented a single cost matrix
for the chemicals, and only a single approximate GED algorithm. While it is unlikely that we will
ever add new matrices or algorithms, we are open to contributions that would do so. If interested
make sure to read through the documentation to understand how to
:ref:`construct new cost matrices <custom_cost_matrix>`.

Setting Up Development Environment
==================================

1. **Fork the repository**

   Start by forking the ChemGED repository on GitHub.

2. **Clone your fork**

   .. code-block:: bash

       git clone https://github.com/YOUR_USERNAME/ChemGED.git
       cd ChemGED

3. **Set up a virtual environment**

   We recommend using Poetry for dependency management:

   .. code-block:: bash

       # Install Poetry if you don't have it
       pip install poetry

       # Create a virtual environment and install dependencies
       poetry install --with dev,test,docs

   Alternatively, you can use pip:

   .. code-block:: bash

       # Create a virtual environment
       python -m venv venv
       source venv/bin/activate  # On Windows: venv\Scripts\activate

       # Install dependencies
       pip install -e ".[dev,test,docs]"

Development Workflow
====================

1. **Create a new branch**

   Always create a new branch for your changes:

   .. code-block:: bash

       git checkout -b feature/your-feature-name

2. **Make your changes**

   Implement your changes, following the coding standards described below.

3. **Run tests**

   Make sure all tests pass:

   .. code-block:: bash

       pytest

   To run tests with coverage:

   .. code-block:: bash

       pytest --cov=chemged

4. **Update documentation**

   If you've added or modified functionality, update the documentation accordingly.

5. **Commit your changes**

   Follow the commit message guidelines below.

6. **Push your changes**

   .. code-block:: bash

       git push origin feature/your-feature-name

7. **Create a pull request**

   Go to the GitHub repository and create a pull request from your branch to the main branch.

Coding Standards
================

We follow standard Python coding conventions:

1. **PEP 8**

   Follow the PEP 8 style guide for Python code. You can use tools like flake8 or black to help with this.

2. **Type Hints**

   Use type hints for function and method signatures:

   .. code-block:: python

       def example_function(param1: str, param2: int) -> bool:
           # ...

3. **Docstrings**

   Use NumPy-style docstrings for all functions, classes, and methods:

   .. code-block:: python

       def example_function(param1, param2):
           """
           Brief description of the function.

           Parameters
           ----------
           param1 : type
               Description of param1
           param2 : type
               Description of param2

           Returns
           -------
           type
               Description of return value
           """
           # ...

4. **Tests**

   Write tests for all new functionality using pytest.
   It is best to make a test case for each new function to maintain high coverage.

ChemGED provides a ``pre-commit`` configuration to help with coding standard checks.
These checks will be enforced prior to merging a pull request, so it is best to check them
before you make a pull results.

To set up pre-commit hooks, run:

.. code-block:: bash

     pre-commit install

This will install the pre-commit hooks defined in the `.pre-commit-config.yaml` file.


Commit Message Guidelines
=========================

We follow the conventional commits specification:

- **feat**: A new feature
- **fix**: A bug fix
- **docs**: Documentation only changes
- **style**: Changes that do not affect the meaning of the code (white-space, formatting, etc)
- **refactor**: A code change that neither fixes a bug nor adds a feature
- **perf**: A code change that improves performance
- **test**: Adding missing tests or correcting existing tests
- **chore**: Changes to the build process or auxiliary tools

Example:

.. code-block:: text

    feat(cost): add new atomic property cost matrix

    This commit adds a new cost matrix implementation that uses atomic properties
    to determine costs for GED calculation.

Pull Request Process
====================

1. Ensure all tests pass.
2. Update the documentation if necessary.
3. The PR should be reviewed by at least one maintainer.
4. Once approved, a maintainer will merge the PR.

Reporting Bugs
==============

If you find a bug, please report it by creating an issue on GitHub. Include the following information:

1. A clear description of the bug
2. Steps to reproduce the bug
3. Expected behavior
4. Actual behavior
5. Environment information (Python version, OS, etc.)
6. If possible, a minimal code example that reproduces the bug

Feature Requests
================

If you have an idea for a new feature, please create an issue on GitHub with the following information:

1. A clear description of the feature
2. The motivation for the feature
3. If possible, a sketch of how the feature might be implemented


Questions?
==========

If you have any questions about contributing, feel free to open an issue on GitHub or reach
out to the maintainers directly.

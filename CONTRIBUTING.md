# Contributing to Kalign

First off, thank you for considering contributing to Kalign! It's people like you that make Kalign such a powerful multiple sequence alignment tool for the bioinformatics community.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [How Can I Contribute?](#how-can-i-contribute)
- [Getting Started](#getting-started)
- [Development Environment Setup](#development-environment-setup)
- [Pull Request Process](#pull-request-process)
- [Issue Guidelines](#issue-guidelines)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Community](#community)

## Code of Conduct

This project and everyone participating in it is governed by our commitment to providing a welcoming and inclusive environment. By participating, you are expected to uphold high standards of respectful communication and collaboration.

**Our Standards:**
- Use welcoming and inclusive language
- Be respectful of differing viewpoints and experiences
- Gracefully accept constructive criticism
- Focus on what is best for the community
- Show empathy towards other community members

## How Can I Contribute?

### Types of Contributions Welcome

We welcome many different types of contributions:

- **Bug reports and fixes** - Help us identify and resolve issues
- **Performance improvements** - Optimizations for speed and memory usage
- **New features** - Alignment algorithms, output formats, or usability enhancements
- **Documentation** - API docs, tutorials, examples, or README improvements
- **Testing** - Writing tests, testing on different platforms, benchmarking
- **Python bindings** - Improvements to the Python package
- **Build system** - CMake, Zig, CI/CD improvements

### What We're NOT Looking For

Please don't use the issue tracker for:
- Support questions (use discussions instead)
- Feature requests that fundamentally change Kalign's core purpose
- Issues related to third-party tools or dependencies

## Getting Started

### Your First Contribution

Unsure where to begin? You can start by looking through these issues:

- **Good first issue** - Issues that are good for newcomers
- **Help wanted** - Issues that need attention from the community

Never made an open source contribution before? Here are some helpful resources:
- [First Timers Only](https://www.firsttimersonly.com/)
- [How to Contribute to an Open Source Project on GitHub](https://opensource.guide/how-to-contribute/)

### Before You Start

For large changes, please open an issue first to discuss what you would like to change. This helps ensure your contribution aligns with the project's goals and avoids duplicate work.

## Development Environment Setup

### Prerequisites

- **C compiler** - GCC, Clang, or MSVC
- **CMake** (3.18 or higher)
- **Git**
- **OpenMP** (optional, for parallelization)

### Optional Dependencies

- **Zig** (for alternative build system)
- **Python** (3.9+ for Python bindings)
- **pybind11** (for Python module development)

### Building from Source

```bash
# Clone your fork
git clone https://github.com/yourusername/kalign.git
cd kalign

# Create build directory
mkdir build && cd build

# Configure and build
cmake ..
make

# Run tests
make test
```

### Python Development

```bash
cd python
pip install -e .
python -c "import kalign; print(kalign.__version__)"
```

## Pull Request Process

1. **Fork** the repository and create your branch from `main`
2. **Make your changes** following our coding standards
3. **Add tests** for any new functionality
4. **Update documentation** if you change APIs or add features
5. **Ensure tests pass** locally before submitting
6. **Submit a pull request** with a clear description

### Pull Request Guidelines

- **One feature per PR** - Keep changes focused and atomic
- **Clear commit messages** - Use descriptive commit messages
- **Link related issues** - Reference any related issue numbers
- **Update CHANGELOG** - Add a brief description of your changes
- **Cross-platform compatibility** - Ensure changes work on Linux, macOS, and Windows

## Issue Guidelines

### Bug Reports

When reporting bugs, please include:

- **Kalign version** - Output of `kalign --version`
- **Operating system** - Including version
- **Build information** - How you installed/built Kalign
- **Input data** - Sample sequences that trigger the bug (if possible)
- **Expected behavior** - What you expected to happen
- **Actual behavior** - What actually happened
- **Error messages** - Complete error output
- **Steps to reproduce** - Minimal steps to trigger the issue

### Feature Requests

When suggesting features:

- **Use case** - Describe the problem you're trying to solve
- **Proposed solution** - How you envision the feature working
- **Alternatives** - Other approaches you've considered
- **Impact** - Who would benefit from this feature

### Performance Issues

For performance problems:

- **Input size** - Number and length of sequences
- **System specs** - CPU, RAM, OS
- **Timing data** - How long operations take
- **Comparison** - Performance with other tools (if applicable)

## Coding Standards

### C Code Style

- **Indentation** - Use spaces, not tabs (consistent with existing code)
- **Naming** - Use descriptive variable and function names
- **Comments** - Document complex algorithms and non-obvious code
- **Memory management** - Always check malloc/free pairs
- **Error handling** - Use consistent error checking patterns

### Python Code Style

- **PEP 8** - Follow Python style guidelines
- **Type hints** - Use type annotations for function signatures
- **Docstrings** - Document all public functions and classes

### Commit Message Format

```
type(scope): brief description

Detailed explanation if needed

Fixes #issue_number
```

Examples:
- `fix(alignment): handle empty sequences correctly`
- `feat(python): add support for custom gap penalties`
- `docs(readme): update installation instructions`

## Testing

### Running Tests

```bash
# C/C++ tests
cd build
make test

# Python tests
cd python
python -m pytest
```

### Test Requirements

- **All tests must pass** before submitting PRs
- **Add tests** for new features and bug fixes
- **Cross-platform testing** - Test on different operating systems when possible
- **Performance tests** - Include benchmarks for performance-critical changes

### Test Data

Use the provided test sequences in `tests/data/` or create minimal test cases. For large datasets, provide instructions for users to download test data separately.

## Community

### Getting Help

- **GitHub Discussions** - For questions and general discussion
- **Issues** - For bug reports and feature requests
- **Email** - Contact the maintainer for security issues

### Recognition

Contributors are recognized in:
- **AUTHORS file** - All contributors are listed
- **Release notes** - Significant contributions are highlighted
- **Paper acknowledgments** - For substantial algorithmic contributions

## License

By contributing to Kalign, you agree that your contributions will be licensed under the GNU General Public License v3.0.

## Questions?

Don't hesitate to ask! The worst thing that can happen is that you'll be politely asked to change something. We appreciate any sort of contribution and don't want a wall of rules to get in the way of that.

Thank you for contributing to Kalign! 🧬
from importlib.util import find_spec


def are_packages_available(packages):
    """Check if packages are available."""
    not_found = [packages[pkg] for pkg in packages if find_spec(pkg) is None]
    not_found = ['scikit-image' if pkg == 'skimage' else pkg for pkg in not_found]
    if len(not_found) != 0:
        msg = f'Packages {", ".join(not_found)} are required for the multimetric submodule. Install before proceeding.'
        raise ImportError(msg)


are_packages_available({'matplotlib': 'matplotlib', 'pandas': 'pandas', 'lmfit': 'lmfit', 'skimage': 'scikit-image'})

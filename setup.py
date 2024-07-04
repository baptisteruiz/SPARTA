from setuptools import setup

setup(
    name='SPARTA',
    version='0.0.1',    
    description='A pipeline for Machine-Learning based analysis of taxonomic profilings of microbial communities',
    url='https://github.com/baptisteruiz/SPARTA.git',
    author='Baptiste Ruiz',
    author_email='baptiste.ruiz@inria.fr',
    license='GNU General Public License v3.0',
    packages=['SPARTA'],
    install_requires=['pandas>=1.4.3',
                        'numpy>=1.21.2',
                        'scikit-learn>=1.1.1',
                        'scipy>=1.7.3',
                        'matplotlib>=3.5.1',
                        'seaborn>=0.12.2',
                        'joblib>=1.1.0',
                        'tqdm',
                        'goatools>=1.2.3',
                        'Biopython>=1.79',
                        'requests>=2.28.1',
                        'kneebow>=1.0.1',
                        'esmecata>=0.2.12',
                        'shap>=0.46.0'                    
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GPL 3.0 License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
)

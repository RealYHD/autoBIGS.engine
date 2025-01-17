pipeline {
    agent {
        kubernetes {
            cloud 'rsys-devel'
            defaultContainer 'miniforge'
            inheritFrom 'miniforge'
        }
    }
    stages {
        stage("install") {
            steps {
                sh 'conda env update -n base -f environment.yml'
            }
        }
        stage("unit tests") {
            steps {
                sh returnStatus: true, script: "pytest --junitxml=test_results.xml --cov=src --cov-report xml:coverage.xml"
                xunit checksName: '', tools: [JUnit(excludesPattern: '', pattern: 'test_results.xml', stopProcessingIfError: true)]
                recordCoverage(tools: [[parser: 'COBERTURA', pattern: 'coverage.xml']])
            }
        }
        stage("build") {
            steps {
                sh "build"
                sh "grayskull pypi dist/*.tar.gz"
                sh "conda-build automlst.engine"
            }
        }
        stage("archive") {
            steps {
                archiveArtifacts artifacts: 'dist/*.tar.gz, dist/*.whl', fingerprint: true, followSymlinks: false, onlyIfSuccessful: true
            }
        }
        stage("publish") {
            parallel {
                stage ("internal") {
                    environment {
                        TOKEN = credentials('git.reslate.systems')
                    }
                    steps {
                        sh returnStatus: true, script: 'python -m twine upload --repository-url https://git.reslate.systems/api/packages/ydeng/pypi -u __token__ -p ${TOKEN} --non-interactive --disable-progress-bar --verbose dist/*'
                    }
                }
                stage ("external") {
                    when {
                        tag '*.*'
                    }
                    environment {
                        PYPI_TOKEN = credentials('pypi.org')
                        CONDA_TOKEN = credentials('anaconda.org')
                    }
                    steps {
                        sh returnStatus: true, script: 'python -m twine upload -u __token__ -p ${TOKEN} --non-interactive --disable-progress-bar --verbose dist/*'
                    }
                }
            }
        }
    }
}
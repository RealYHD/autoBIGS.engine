pipeline {
    agent {
        kubernetes {
            cloud 'rsys-devel'
            defaultContainer 'miniforge3'
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
            environment {
                GIT_AUTHOR_NAME = "Harrison Deng"
                GIT_AUTHOR_EMAIL = "yunyangdeng@outlook.com"
            }
            steps {
                sh "python -m build"
                sh "grayskull pypi dist/*.tar.gz"
                sh "conda-build autobigsst.engine --output-folder conda-bld"
            }
        }
        stage("archive") {
            steps {
                archiveArtifacts artifacts: 'dist/*.tar.gz, dist/*.whl conda-bld/**/*.conda', fingerprint: true, followSymlinks: false, onlyIfSuccessful: true
            }
        }
        stage("publish") {
            parallel {
                stage ("internal") {
                    environment {
                        CREDS = credentials('username-password-rs-git')
                    }
                    steps {
                        sh returnStatus: true, script: 'python -m twine upload --repository-url https://git.reslate.systems/api/packages/ydeng/pypi -u ${CREDS_USR} -p ${CREDS__PSW} --non-interactive --disable-progress-bar --verbose dist/*'
                    }
                }
                stage ("external") {
                    when {
                        tag '*.*'
                    }
                    environment {
                        PYPI_TOKEN = credentials('pypi.org')
                    }
                    steps {
                        sh returnStatus: true, script: 'python -m twine upload -u __token__ -p ${TOKEN} --non-interactive --disable-progress-bar --verbose dist/*'
                    }
                }
            }
        }
    }
}
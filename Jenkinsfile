pipeline {
    agent {
        kubernetes {
            cloud 'rsys-devel'
            defaultContainer 'pip'
            inheritFrom 'pip'
        }
    }
    stages {
        stage("install") {
            steps {
                sh 'python -m pip install -r requirements.txt'
            }
        }
        stage("unit tests") {
            steps {
                sh returnStatus: true, script: "python -m pytest --junitxml=test_results.xml --cov=src --cov-report xml:coverage.xml"
                xunit checksName: '', tools: [JUnit(excludesPattern: '', pattern: 'test_results.xml', stopProcessingIfError: true)]
                recordCoverage(tools: [[parser: 'COBERTURA', pattern: 'coverage.xml']])
            }
        }
        stage("build") {
            steps {
                sh "python -m build"
            }
        }
        stage("archive") {
            steps {
                archiveArtifacts artifacts: 'dist/*.tar.gz, dist/*.whl', fingerprint: true, followSymlinks: false, onlyIfSuccessful: true
            }
        }
        stage("publish") {
            parallel {
                stage ("git.reslate.systems") {
                    environment {
                        TOKEN = credentials('git.reslate.systems')
                    }
                    steps {
                        sh returnStatus: true, script: 'python -m twine upload --repository-url https://git.reslate.systems/api/packages/ydeng/pypi -u __token__ -p ${TOKEN} --non-interactive --disable-progress-bar --verbose dist/*'
                    }
                }
                stage ("pypi.org") {
                    when {
                        tag '*.*'
                    }
                    environment {
                        TOKEN = credentials('pypi.org')
                    }
                    steps {
                        sh returnStatus: true, script: 'python -m twine upload -u __token__ -p ${TOKEN} --non-interactive --disable-progress-bar --verbose dist/*'
                    }
                }
            }
        }
    }
}
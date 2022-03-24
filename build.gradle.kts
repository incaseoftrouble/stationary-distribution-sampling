plugins {
    java
    application
    idea
}

group = "mc"
version = "0.1"

base {
    archivesName.set("stationary-distribution-sampling")
}

java {
    sourceCompatibility = JavaVersion.VERSION_17
    targetCompatibility = JavaVersion.VERSION_17

    withSourcesJar()
    withJavadocJar()
}

var defaultEncoding = "UTF-8"
tasks.withType<JavaCompile> { options.encoding = defaultEncoding }
tasks.withType<Javadoc> { options.encoding = defaultEncoding }
tasks.withType<Test> { systemProperty("file.encoding", "UTF-8") }

tasks.test {
    useJUnitPlatform()
}

idea {
    module {
        isDownloadJavadoc = true
        isDownloadSources = true
    }
}

repositories {
    mavenCentral()
}

val modelsProject = project(":lib:models")
dependencies {
    implementation(modelsProject)
    implementation(files("lib/Jeigen-onefat.jar"))
    implementation("net.java.dev.jna", "jna", "4.1.0")
    implementation("info.picocli", "picocli", "4.6.2")

    testImplementation("org.junit.jupiter", "junit-jupiter-api", "5.8.2")
    testRuntimeOnly("org.junit.jupiter", "junit-jupiter-engine", "5.8.2")
}

application {
    mainClass.set("stationary.Main")
    // applicationDefaultJvmArgs = listOf("-Djava.library.path=build/prism/natives")
}

val startScripts = tasks.getByName("startScripts", CreateStartScripts::class)

tasks.register("extractTemplate") {
    file(rootProject.buildDir).mkdirs()
    file("${rootProject.buildDir}/template-unix.txt").writeText(
        (startScripts.unixStartScriptGenerator as TemplateBasedScriptGenerator).template.asString()
    )
}

startScripts.doFirst {
    (startScripts.unixStartScriptGenerator as TemplateBasedScriptGenerator).template =
        project.resources.text.fromFile("${modelsProject.projectDir}/config/template-unix.txt")
}
startScripts.doLast { startScripts.windowsScript.delete() }

distributions {
    main {
        contents {
            modelsProject.afterEvaluate {
                from(modelsProject.tasks.getByName("prismNatives")) {
                    into("lib")
                }
            }
        }
    }
}

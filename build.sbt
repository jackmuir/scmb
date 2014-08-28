name := "cmb_hmc"

version := "0.1"

libraryDependencies  ++= Seq(
            // other dependencies here
            "org.scalanlp" % "breeze_2.10" % "0.9",
            // native libraries are not included by default. add this if you want them (as of 0.7)
            // native libraries greatly improve performance, but increase jar sizes.
            "org.scalanlp" % "breeze-natives_2.10" % "0.9"
)

resolvers ++= Seq(
            // other resolvers here
            "Sonatype Releases" at "https://oss.sonatype.org/content/repositories/releases/"
)


// Scala 2.9.2 is still supported for 0.2.1, but is dropped afterwards.
// Don't use an earlier version of 2.10, you will probably get weird compiler crashes.
scalaVersion := "2.11.2"
